#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <fstream>

#include <filesystem>
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "overlap_source.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

class RamOverlapSource : public OverlapSource {
private:
    std::vector<std::vector<biosoup::Overlap>> overlaps;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    std::string reads;
    unsigned int kmer = 15; // 30
    unsigned int w = 5; // 15
    unsigned int b = 500;
    unsigned int chain = 4;
    unsigned int matches = 100;
    unsigned int gap = 10000;
    double f = 0.001;
    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> create_parser(std::string& path){
        auto is_suffix = [] (const std::string& s, const std::string& suff) {
            return s.size() >= suff.size() && s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
        };

        if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
            is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz")   ||
            is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz"))  {
            try {
                return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);

            } catch (const std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
            }
        }
        return nullptr;
    }

    void load_sequences(){
        auto parser = create_parser(reads);
        while (true){
            std::vector<std::unique_ptr<biosoup::NucleicAcid>> buffer;
            try{
                buffer = parser->Parse(1U <<30);
            }catch (std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
                exit(1);
            }

            if(buffer.empty()){
                break;
            }
            sequences.reserve(sequences.size() + buffer.size());
            for(const auto& sequence: buffer){
                sequences.push_back(std::make_unique<biosoup::NucleicAcid>(*sequence));
            }
        }

    }

    void generate(){
        auto threads = std::make_shared<thread_pool::ThreadPool>(64);
        ram::MinimizerEngine minimizer_engine{
                threads,
                kmer,
                w,
                b,
                chain,
                matches,
                gap
        };

        std::size_t bytes = 0;
        for(std::uint32_t  i = 0, j=0; i < sequences.size(); i++){
            bytes += sequences[i]->inflated_len;
            if (i != sequences.size() - 1 && bytes < (1ULL << 32)) {
                continue;
            }
            bytes = 0;
            minimizer_engine.Minimize(
                    sequences.begin() + j,
                    sequences.begin() + i + 1,
                    true);
            minimizer_engine.Filter(f);

            std::vector<std::uint32_t> num_overlaps(overlaps.size());
            for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
                num_overlaps[k] = overlaps[k].size();
            }

            std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;

            auto edlib_wrapper = [&](
                    std::uint32_t i,
                    const biosoup::Overlap &it,
                    const std::string &lhs,
                    const std::string &rhs) -> std::string {
                std::string cigar_alignment = "";
                EdlibAlignResult result = edlibAlign(
                        lhs.c_str(), lhs.size(),
                        rhs.c_str(), rhs.size(),
                        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0)); // align lhs and rhs
                if (result.status == EDLIB_STATUS_OK) {
                    cigar_alignment = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
                    edlibFreeAlignResult(result);
                    return cigar_alignment;
                } else {
                    edlibFreeAlignResult(result);
                    std::string cigar_alignment = "";
                    return cigar_alignment;
                }
            };
            auto cigar_alignment_reverse = [](const std::string &s) -> std::string {
                std::string rs;
                if (s.empty()) {
                    return rs;
                }
                for (char c: s) {
                    switch (c) {
                        case 'I':
                            rs += 'D';
                            break;
                        case 'D':
                            rs += 'I';
                            break;
                        default:
                            rs += c;
                            break;
                    }
                }
                return rs;
            };

            auto cigar_overlap_reverse = [&cigar_alignment_reverse](const biosoup::Overlap &o) -> biosoup::Overlap {
                return biosoup::Overlap(
                        o.rhs_id, o.rhs_begin, o.rhs_end,
                        o.lhs_id, o.lhs_begin, o.lhs_end,
                        o.score, cigar_alignment_reverse(o.alignment),
                        o.strand);
            };

            auto overlap_length = [](const biosoup::Overlap &o) -> std::uint32_t {
                return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
            };

            for (std::uint32_t k = 0; k < i + 1; ++k) {
                thread_futures.emplace_back(threads->Submit(
                        [&](std::uint32_t i) -> std::vector<biosoup::Overlap> { // map sequences and fill out the potential snp list
                            std::vector<biosoup::Overlap> ovlps = minimizer_engine.Map(sequences[i], true, true,
                                                                                       true, false);
                            if (!ovlps.empty()) {
                                std::vector<biosoup::Overlap> ovlps_final;

                                std::sort(ovlps.begin(), ovlps.end(),
                                          [&](const biosoup::Overlap &lhs,
                                              const biosoup::Overlap &rhs) -> bool {
                                              return overlap_length(lhs) > overlap_length(rhs);
                                          });

                                std::vector<biosoup::Overlap> tmp;
                                tmp.insert(tmp.end(), ovlps.begin(),
                                           ovlps.begin() + (ovlps.size() > 24 ? 24 : ovlps.size()));  // NOLINT
                                tmp.swap(ovlps);

                                for (const auto &ovlp: ovlps) {
                                    auto lhs = sequences[i]->InflateData(ovlp.lhs_begin,
                                                                         ovlp.lhs_end - ovlp.lhs_begin);
                                    biosoup::NucleicAcid rhs_{"",
                                                              sequences[ovlp.rhs_id]->InflateData(ovlp.rhs_begin,
                                                                                                  ovlp.rhs_end -
                                                                                                  ovlp.rhs_begin)};

                                    if (!ovlp.strand) rhs_.ReverseAndComplement();

                                    auto rhs = rhs_.InflateData();

                                    std::string tmp = edlib_wrapper(i, ovlp, lhs, rhs);

                                    biosoup::Overlap ovlp_tmp{ovlp.lhs_id, ovlp.lhs_begin, ovlp.lhs_end,
                                                              ovlp.rhs_id, ovlp.rhs_begin,
                                                              ovlp.rhs_end, ovlp.score, tmp, ovlp.strand};
                                    ovlps_final.emplace_back(ovlp_tmp);

                                };

                                return ovlps_final;
                            }
                            return ovlps;
                        },
                        k));

                bytes += sequences[k]->inflated_len;
                if (k != i && bytes < (1U << 30)) {
                    continue;
                }
                bytes = 0;

                for (auto &it: thread_futures) {
                    std::vector<biosoup::Overlap> overlaps_tmp;
                    for (const auto &jt: it.get()) {
                        overlaps_tmp.emplace_back(jt);
                        //overlaps[jt.lhs_id].emplace_back(jt);
                        //overlaps[jt.rhs_id].emplace_back(cigar_overlap_reverse(jt));

                    }
                    std::sort(overlaps_tmp.begin(), overlaps_tmp.end(), [](auto  o1, auto o2) {
                        return o1.score > o2.score;
                    });
                    for(size_t v = 0; v < std::min(30UL,overlaps_tmp.size()); v++){
                        overlaps[overlaps_tmp[v].lhs_id].emplace_back(overlaps_tmp[v]);
                        overlaps[overlaps_tmp[v].rhs_id].emplace_back(cigar_overlap_reverse(overlaps_tmp[v]));

                    }
                }
                thread_futures.clear();

                std::vector<std::future<void>> void_futures;

                for (const auto &it: void_futures) {
                    it.wait();
                }
            }
            j = i + 1;
        }
    }

public:
    RamOverlapSource(const std::string& reads){
        this->reads = reads;
        load_sequences();
        this->overlaps.resize(sequences.size());
    }
    std::vector<std::vector<biosoup::Overlap>>* get_overlaps() override {
        if(overlaps.empty()) generate();
        return &overlaps;
    }
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> * get_sequences() override{
        return &sequences;
    }


};

extern "C" OverlapSource* __attribute__((visibility("default"))) create(std::string& args){
    return (OverlapSource*) new RamOverlapSource(args);
}