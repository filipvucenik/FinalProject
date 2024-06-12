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
#include "overlap_sources/overlap_source_factory.hpp"
#include "thread_pool/thread_pool.hpp"
#include "edlib.h"



std::uint8_t use_frequencies = 1;
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


struct snp_position {
    size_t sequence_ind;
    size_t position;
    base_pile pile;
};
struct informative_position{
    std::uint32_t informative;
    std::uint32_t informative_matches;
    std::uint32_t informative_mismatches;
};

void printPaf(std::vector<std::vector<biosoup::Overlap>>* overlaps, std::vector<std::unique_ptr<biosoup::NucleicAcid>>* sequences, std::string& out_file, std::vector<std::vector<informative_position>> &informative_positions){
    std::ofstream os(out_file);
    std::uint32_t i = 0;
    for (const auto &it: *overlaps) {
        std::uint32_t j = 0;
        for (const auto &jt: it) {
            std::string positions = "";
            bool haplotypes = false;
            if(informative_positions[i][j].informative != 0 && informative_positions[i][j].informative_matches){
                haplotypes = (double)informative_positions[i][j].informative_matches / informative_positions[i][j].informative < 0.8;
            }
            positions = std::to_string(informative_positions[i][j].informative_matches) + "/" + std::to_string(informative_positions[i][j].informative);
            os << (*sequences)[jt.lhs_id]->name
               << "\t" << (*sequences)[jt.lhs_id]->inflated_len  // length
               << "\t" << jt.lhs_begin
               << "\t" << jt.lhs_end
               << "\t" << (jt.strand ? "+" : "-")
               << "\t" << (*sequences)[jt.rhs_id]->name
               << "\t" << (*sequences)[jt.rhs_id]->inflated_len  // length
               << "\t" << jt.rhs_begin
               << "\t" << jt.rhs_end
               << "\t" << 255 // residue matches
               << "\t" << 255 // alignment block length
               << "\t" << jt.alignment
               << "\t" << haplotypes
               << "\t" << positions
               << std::endl;
            j++;
        }
        i++;
    }
}




int main (int argc, char* argv[]) {

    std::filesystem::path currentDir = std::filesystem::current_path();
    if (argc < 3) {
        std::cout << usage;
        exit(1);
    }
    std::string  lib = argv[1];
    std::string arg = argv[2];
    if(argc > 3){
        for(int i = 3 ; i < argc; i++){
            arg += " ";
            arg += argv[i];
        }
    }

    OverlapSource* os = create_overlap_source(lib, arg);
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>* sequences = os->get_sequences();
    std::cout<<"sequences loaded sequences"<<(*sequences).size()<<std::endl;
    std::vector<std::vector<biosoup::Overlap>>* overlaps = os->get_overlaps();
    std::vector<std::unordered_set<std::uint32_t>> annotations((*sequences).size());
    std::vector<std::vector<snp_position>> snp_positions(50);

    std::vector<std::vector<informative_position>> informative_positions((*overlaps).size());



    for(std::size_t i = 0; i < informative_positions.size(); i++){
        for(size_t j = 0; j < (*overlaps)[i].size(); j++){
            informative_positions[i].emplace_back();
        }
    }
    std::cout<<"setup finished overlaps "<<(*overlaps).size()<<std::endl;

    auto cigar_to_edlib_alignment = [](const std::string &s) -> std::string {
        std::string rs = "";
        std::uint64_t pos = 0;
        std::uint64_t start_pos = 0;
        std::uint64_t total_num = 0;

        if (s.empty()) {
            return rs;
        }
        for (size_t i = 0; i < s.length(); i++) {
            if (std::isdigit(s[i])) {
                if (pos == 0) {
                    start_pos = i;
                }
                ++pos;
            } else {
                total_num = 0;
                for (size_t j = start_pos; j < start_pos + pos; j++) {
                    total_num += (s[j] - '0') * std::pow(10, (start_pos + pos) - j - 1);
                }
                pos = 0;
                switch (s[i]) {
                    case '=':
                        for (size_t j = 0; j < total_num; j++) {
                            rs += '\000';
                        };
                        break;
                    case 'X':
                        for (size_t j = 0; j < total_num; j++) {
                            rs += '\003';
                        };
                        break;
                    case 'I':
                        for (size_t j = 0; j < total_num; j++) {
                            rs += '\001';
                        };
                        break;
                    case 'D':
                        for (size_t j = 0; j < total_num; j++) {
                            rs += '\002';
                        };
                        break;
                    case 'M':
                        for(size_t j = 0; j < total_num; j++){
                            rs += '\000';
                        }
                    default:
                        //rs += '\000';
                        break;
                }
            }
        }
        return rs;
    };
    auto call_snps = [&](std::uint32_t i, std::vector<biosoup::Overlap>& ovlps_final) -> void {
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
                            if (lhs_pos < base_pile_tmp.size()) {
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

        std::uint32_t o = 0;
        for(auto &ovlp: ovlps_final){
            std::uint32_t counter_informative = 0;
            std::uint32_t counter_informative_matches = 0;
            std::uint32_t counter_informative_mismatches = 0;
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
                        case 0:{
                            std::uint32_t most_occuring = std::max(base_pile_tmp[lhs_pos].a, std::max(base_pile_tmp[lhs_pos].c, std::max(base_pile_tmp[lhs_pos].g, base_pile_tmp[lhs_pos].t)));
                            std::uint32_t pile_size = base_pile_tmp[lhs_pos].a + base_pile_tmp[lhs_pos].c + base_pile_tmp[lhs_pos].g + base_pile_tmp[lhs_pos].t;
                            if((double) most_occuring / pile_size < 0.8){
                                counter_informative++;
                                counter_informative_matches++;
                            }

                            ++lhs_pos;
                            ++rhs_pos;
                            break;
                        }
                        case 3: {
                            std::uint32_t most_occuring = std::max(base_pile_tmp[lhs_pos].a, std::max(base_pile_tmp[lhs_pos].c, std::max(base_pile_tmp[lhs_pos].g, base_pile_tmp[lhs_pos].t)));
                            std::uint32_t pile_size = base_pile_tmp[lhs_pos].a + base_pile_tmp[lhs_pos].c + base_pile_tmp[lhs_pos].g + base_pile_tmp[lhs_pos].t;
                            if((double) most_occuring / pile_size < 0.8){
                                counter_informative++;
                                counter_informative_mismatches++;
                            }
                            ++lhs_pos;
                            ++rhs_pos;
                            break;
                        }
                        case 1: {
                            ++base_pile_tmp[lhs_pos].i;
                            ++lhs_pos;
                            break;
                        }
                        case 2: {
                            if (lhs_pos < base_pile_tmp.size()) {
                                ++base_pile_tmp[lhs_pos].d;
                            }
                            ++rhs_pos;
                            break;
                        }
                        default:
                            break;
                    }
                }
                informative_positions[i][o].informative = counter_informative;
                informative_positions[i][o].informative_matches = counter_informative_matches;
                informative_positions[i][o].informative_mismatches = counter_informative_mismatches;
            }
            o++;
        }

        cov.reserve(base_pile_tmp.size());
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
            if(i < 50){
                snp_positions[i].push_back({i, j, jt});
            }

            double sum = std::accumulate(counts.begin(), counts.end(), 0);

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
    auto threads = std::make_shared<thread_pool::ThreadPool>(64);
    std::vector<std::future<void>> futures;

    for(size_t l = 0; l < (*sequences).size(); l++){
        call_snps(l, (*overlaps)[l]);
        /*
        futures.emplace_back(threads->Submit(
                [&](size_t l){
                    call_snps(l, (*overlaps)[l]);
                },l
                ));
        */
    }

    for(auto &future: futures){
        future.wait();
    }

    std::string out_file = "hifiasm_overlapsPA_nanosim.paf";
    printPaf(overlaps, sequences, out_file, informative_positions);


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

    /*
    std::ofstream pile_file;
    pile_file.open("pile.txt");
    */


    for(auto &it: snp_positions){
        if(it.empty()){
            continue;
        }
        size_t ind = it[0].sequence_ind;
        std::string seq = (*sequences)[ind]->InflateData();
        /*
        pile_file <<"sequence: "<< ind << std::endl;
       // pile_file<<seq<<std::endl;
        for(auto jt: it){
            pile_file<<seq[jt.position]<<"|a:"<<jt.pile.a<<"c:"<<jt.pile.c<<"g:"<<jt.pile.g<<"t:"<<jt.pile.t<<"i:"<<jt.pile.i<<"d:"<<jt.pile.d<<std::endl;
        }
         */
    }
    return 0;

}