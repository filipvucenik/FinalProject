#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <fstream>

#include <filesystem>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "overlap_sources/overlap_source_factory.hpp"



std::uint8_t use_frequencies = 0;
std::uint8_t variant_call_th = 3;
double freq_low_th = 0.333;
double freq_high_th = 0.667;
std::uint8_t print_snp_data = 1;


std::string usage = "Usage: snp <overlap_source> <overlap_source_args>\n";



struct base_pile {
    std::uint32_t a;
    std::uint32_t c;
    std::uint32_t g;
    std::uint32_t t;
    std::uint32_t i;
    std::uint32_t d;
};




int main (int argc, char* argv[]) {
    std::filesystem::path currentDir = std::filesystem::current_path();
    if (argc != 2) {
        std::cout << usage;
        exit(1);
    }
    std::string reads_file = argv[1];

    OverlapSource* os = create_overlap_source("ram", reads_file);
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>* sequences = os->get_sequences();
    std::vector<std::vector<biosoup::Overlap>>* overlaps = os->get_overlaps();
    std::vector<std::unordered_set<std::uint32_t>> annotations((*sequences).size());

    auto overlap_length = [](const biosoup::Overlap &o) -> std::uint32_t {
        return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
    };

    auto edlib_alignment_reverse = [](const std::string &s) -> std::string {
        std::string rs;
        if (s.empty()) {
            return rs;
        }
        for (char c: s) {
            switch (c) {
                case '2':
                    rs += '\001';
                    break;
                case '1':
                    rs += '\002';
                    break;
                default:
                    rs += c;
                    break;
            }
        }
        return rs;
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


    auto cigar_to_edlib_alignment = [](const std::string &s) -> std::string {
        std::string rs = "";
        std::uint64_t pos = 0;
        std::uint64_t start_pos = 0;
        std::uint64_t total_num = 0;

        if (s.empty()) {
            return rs;
        }
        for (int i = 0; i < s.length(); i++) {
            if (std::isdigit(s[i])) {
                if (pos == 0) {
                    start_pos = i;
                }
                ++pos;
            } else {
                total_num = 0;
                for (int j = start_pos; j < start_pos + pos; j++) {
                    total_num += (s[j] - '0') * std::pow(10, (start_pos + pos) - j - 1);
                }
                pos = 0;
                switch (s[i]) {
                    case '=':
                        for (int j = 0; j < total_num; j++) {
                            rs += '\000';
                        };
                        break;
                    case 'X':
                        for (int j = 0; j < total_num; j++) {
                            rs += '\003';
                        };
                        break;
                    case 'I':
                        for (int j = 0; j < total_num; j++) {
                            rs += '\001';
                        };
                        break;
                    case 'D':
                        for (int j = 0; j < total_num; j++) {
                            rs += '\002';
                        };
                        break;
                    default:
                        //rs += '\000';
                        break;
                }
            }
        }
        return rs;
    };

    auto overlap_reverse = [&edlib_alignment_reverse](const biosoup::Overlap &o) -> biosoup::Overlap {
        return biosoup::Overlap(
                o.rhs_id, o.rhs_begin, o.rhs_end,
                o.lhs_id, o.lhs_begin, o.lhs_end,
                o.score, edlib_alignment_reverse(o.alignment),
                o.strand);
    };

    auto cigar_overlap_reverse = [&cigar_alignment_reverse](const biosoup::Overlap &o) -> biosoup::Overlap {
        return biosoup::Overlap(
                o.rhs_id, o.rhs_begin, o.rhs_end,
                o.lhs_id, o.lhs_begin, o.lhs_end,
                o.score, cigar_alignment_reverse(o.alignment),
                o.strand);
    };


    auto call_snps = [&](std::uint32_t i, std::vector<biosoup::Overlap> ovlps_final) -> void {
        std::uint32_t seq_inflated_len = (*sequences)[i]->inflated_len;
        std::vector<base_pile> base_pile_tmp(seq_inflated_len);
        std::vector<std::uint32_t> cov;

        for (auto &ovlp: ovlps_final) {
            if (!(ovlp.alignment.empty())) {
                std::uint32_t lhs_pos = ovlp.lhs_begin;
                std::uint32_t rhs_pos = 0;
                biosoup::NucleicAcid rhs{"",
                                         (*sequences)[ovlp.rhs_id]->InflateData(ovlp.rhs_begin,
                                                                             ovlp.rhs_end - ovlp.rhs_begin)};
                if (!ovlp.strand) rhs.ReverseAndComplement();
                std::string rhs_tmp = rhs.InflateData();

                std::string edlib_alignment = cigar_to_edlib_alignment(ovlp.alignment);

                for (auto &edlib_align: edlib_alignment) {
                    switch (edlib_align) {
                        case 0:
                        case 3: {
                            switch (rhs_tmp[rhs_pos]) {
                                case 'A':
                                    ++base_pile_tmp[lhs_pos].a;
                                    break;
                                case 'C':
                                    ++base_pile_tmp[lhs_pos].c;
                                    break;
                                case 'G':
                                    ++base_pile_tmp[lhs_pos].g;
                                    break;
                                case 'T':
                                    ++base_pile_tmp[lhs_pos].t;
                                    break;
                                default:
                                    break; // if they align
                            }
                            ++lhs_pos;
                            ++rhs_pos;
                            break;
                        }
                        case 1: {
                            ++base_pile_tmp[lhs_pos].i;
                            ++lhs_pos;
                            break; // insertion on the left hand side
                        }
                        case 2: {
                            if (!(lhs_pos >= base_pile_tmp.size())) {
                                ++base_pile_tmp[lhs_pos].d;
                            }
                            ++rhs_pos;
                            break; // deletion on the left hand side
                        }
                        default:
                            break;
                    }
                }
            }
        }

        for (const auto &jt: base_pile_tmp) {
            //cov.emplace_back(jt.a + jt.c + jt.g + jt.t);
            cov.emplace_back(jt.a + jt.c + jt.g + jt.t + jt.d + jt.i);
        }

        std::nth_element(cov.begin(), cov.begin() + cov.size() / 2, cov.end());
        double m = cov[cov.size() / 2] * 2. / 3.;

        std::size_t j = 0;
        for (const auto &jt: base_pile_tmp) {
            std::vector<double> counts = {
                    static_cast<double>(jt.a),
                    static_cast<double>(jt.c),
                    static_cast<double>(jt.g),
                    static_cast<double>(jt.t),
                    static_cast<double>(jt.d),
                    static_cast<double>(jt.i)
            };

            double sum = std::accumulate(counts.begin(), counts.end(), 0);

/*
            if(sum > 0)
                std::cout<<"sum="<<sum<<std::endl;
*/

            if (use_frequencies) {
                for (auto &kt: counts) {
                    kt /= sum;
                };
            };

            if (sum > m) {
                std::size_t variants = 0;
                for (const auto &it: counts) {
                    if (use_frequencies) {
                        if (freq_low_th < it && it < freq_high_th) {
                            ++variants;
                        }
                    } else {
                        if (it > variant_call_th) {
                            ++variants;
                        }
                    }
                }
                if (variants > 1) annotations[i].emplace(j);

            };
            ++j;
        };
    };


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



    for(size_t l = 0; l < (*sequences).size(); l++){
        call_snps(l, (*overlaps)[l]);
    }


    std::cout<<annotations.size()<<std::endl;

    std::ofstream outdata;
    outdata.open("snp_annotations.anno");

    for (std::uint32_t i = 0; i < annotations.size(); ++i) {
        if (annotations[i].empty()) {
            continue;
        }
        outdata << i;
        for (const auto &jt: annotations[i]) {
            outdata << " " << jt;
        }
        outdata << std::endl;
    }



    return 0;

}