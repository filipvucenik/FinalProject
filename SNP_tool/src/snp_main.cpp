#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <fstream>
#include <hash_map>
#include <filesystem>
#include <set>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "overlap_sources/overlap_source_factory.hpp"
#include "thread_pool/thread_pool.hpp"
#include "edlib.h"

#include "overlap_sources/overlap_source.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

enum class filtering
{
    NONE,
    MEAN,
    MEDIAN,
    DOUBLE_PERFECT,
    DOUBLE_TRESHOLD
};

filtering filter = filtering::DOUBLE_TRESHOLD;

bool use_frequencies = true;
std::uint8_t variant_call_th = 1;
double freq_low_th = 0.333;
double freq_high_th = 0.667;
bool indels = false;
size_t cutofparam = 2.5;
bool print_paf = false;
double error_rate = 0.003;
size_t no_snps_thr = 5;
bool hifiasm_output = true;

std::string usage = "Usage: snp <overlap_source> <overlap_source_args>\n";

struct base_pile
{
    std::uint32_t a;
    std::uint32_t c;
    std::uint32_t g;
    std::uint32_t t;
    std::uint32_t i;
    std::uint32_t d;
};

struct snp_position
{
    size_t sequence_ind;
    size_t position;
    base_pile pile;
};
struct informative_position
{
    std::uint32_t informative;
    std::uint32_t informative_matches;
    std::uint32_t informative_mismatches;
};

void printPaf(std::vector<std::vector<biosoup::Overlap>> *overlaps, std::vector<std::unique_ptr<biosoup::NucleicAcid>> *sequences, std::string &out_file)
{
    std::ofstream os(out_file);
    std::uint32_t i = 0;
    for (const auto &it : *overlaps)
    {
        std::uint32_t j = 0;
        for (const auto &jt : it)
        {
            os << (*sequences)[jt.lhs_id]->name
               << "\t" << (*sequences)[jt.lhs_id]->inflated_len // length
               << "\t" << jt.lhs_begin
               << "\t" << jt.lhs_end
               << "\t" << (jt.strand ? "+" : "-")
               << "\t" << (*sequences)[jt.rhs_id]->name
               << "\t" << (*sequences)[jt.rhs_id]->inflated_len // length
               << "\t" << jt.rhs_begin
               << "\t" << jt.rhs_end
               << "\t" << jt.score
               << "\t" << jt.alignment
               << std::endl;
            j++;
        }
        i++;
    }
}

void printPaf(std::vector<std::vector<biosoup::Overlap>> *overlaps, std::unordered_map<size_t, std::vector<biosoup::Overlap>::iterator> &end_filter, std::vector<std::unique_ptr<biosoup::NucleicAcid>> *sequences, std::string &out_file)
{
    std::ofstream os(out_file);
    std::uint32_t i = 0;
    for (const auto &it : *overlaps)
    {
        if (it.empty())
            continue;
        std::uint32_t j = 0;
        for (auto jt = it.begin(); jt != end_filter[it[0].lhs_id]; jt++)
        {
            os << (*sequences)[(*jt).lhs_id]->name
               << "\t" << (*sequences)[(*jt).lhs_id]->inflated_len // length
               << "\t" << (*jt).lhs_begin
               << "\t" << (*jt).lhs_end
               << "\t" << ((*jt).strand ? "+" : "-")
               << "\t" << (*sequences)[(*jt).rhs_id]->name
               << "\t" << (*sequences)[(*jt).rhs_id]->inflated_len // length
               << "\t" << (*jt).rhs_begin
               << "\t" << (*jt).rhs_end
               << "\t" << (*jt).score
               << "\t" << (*jt).alignment
               << std::endl;
            j++;
        }
        i++;
    }
}

int main(int argc, char *argv[])
{
    std::filesystem::path currentDir = std::filesystem::current_path();
    if (argc < 4)
    {
        std::cout << usage;
        exit(1);
    }
    std::string lib = argv[1];
    std::string arg = argv[2];
    std::string sample = argv[argc - 1];
    if (argc > 3)
    {
        for (int i = 3; i < argc - 1; i++)
        {
            arg += " ";
            arg += argv[i];
        }
    }
    std::cout << "currentDir: " << currentDir << std::endl;
    std::cout << "arg: " << arg << std::endl;
    std::cout << "lib: " << lib << std::endl;
    std::cout << "sample: " << sample << std::endl;

    // std::string arg = "../../../../reads/E-coli_reads_15c.fasta ../../../../hifiasmOverlaps/hifiasm_overlaps_15c.paf";
    // std::string lib = "pafWFA2";
    // std::string sample = "E-Coli";

    // std::string arg = "../../../../reads/chr19_perfect.fasta ../../../../hifiasmOverlaps/chr19_ovlp_perfect.paf";
    // std::string lib = "pafEDL";
    // std::string sample = "chr19_perfect";

    // debugging

    // std::string arg = "/home/stufilip/scratch/FinalProject/filtered/reads/chr19_e_EDL.fasta /home/stufilip/scratch/FinalProject/filtered/overlaps/chr19_e_EDL.paf";
    // std::string lib = "pafWFA2";
    // std::string sample = "chr19_perfect_filtered";

    // std::string arg = "../../../../reads/chr19_perfect.fasta /home/stufilip/scratch/FinalProject/PA/chr19_p/hifiasm_overlapsPA_chr19_perfect_K62_pafKSW2.paf";
    // std::string arg = "../../../../reads/chr19_perfect.fasta /home/stufilip/scratch/FinalProject/PA/chr19_p/hifiasm_overlapsPA_chr19_perfect_pafEDL.paf";
    // std::string arg = "../../../../reads/badread_50x_corrected_fixed.fasta /home/stufilip/scratch/FinalProject/SNP_tool/build/bin/hifiasm_overlapsPA_badread_50X_pafEDL.paf";
    // std::string arg = "/home/stufilip/scratch/FinalProject/filtered/reads/aad_b40_debug.fasta /home/stufilip/scratch/FinalProject/filtered/overlaps/aad_b40_debug.paf";
    // std::string lib = "pafPA";
    // std::string sample = "badread_50X_EDL_dbg";

    std::unique_ptr<OverlapSource>
        os(create_overlap_source(lib, arg));
    std::cout << "using library: " << lib << std::endl;

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> *sequences = os->get_sequences();
    std::cout << "sequences loaded sequences" << (*sequences).size() << std::endl;
    std::unique_ptr<std::vector<std::vector<biosoup::Overlap>>> overlaps = os->get_overlaps();

    std::string out_file = "hifiasm_overlapsPA_" + sample + "_" + lib + ".paf";
    std::string out_file_filtered = "hifiasm_overlapsPA_filtered_coverage_" + sample + "_" + lib + ".paf";
    if (print_paf)
    {
        std::cout
            << "Writnig to file: " << out_file << std::endl;
        printPaf(overlaps.get(), sequences, out_file);
    }
    std::vector<std::pair<std::vector<std::pair<std::uint32_t, char>>, double>> annotations((*sequences).size());
    std::vector<std::vector<std::tuple<short, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>>> classified_ovlps((*overlaps).size());

    std::size_t ovlp_sum = std::accumulate((*overlaps).begin(), (*overlaps).end(), 0, [](std::size_t sum, const std::vector<biosoup::Overlap> &ovlps)
                                           { return sum + ovlps.size(); });

    double average_overlaps = ovlp_sum / (double)sequences->size();
    std::cout << "estimated coverage: " << average_overlaps << std::endl;

    std::vector<std::size_t> ovlps_per_read = std::vector<std::size_t>(sequences->size());
    for (std::size_t i = 0; i < (*overlaps).size(); i++)
    {
        ovlps_per_read[i] = (*overlaps)[i].size();
    }
    std::sort(ovlps_per_read.begin(), ovlps_per_read.end());
    std::size_t median_overlaps = ovlps_per_read[ovlps_per_read.size() / 2];

    std::cout << "median overlaps: " << median_overlaps << std::endl;

    std::unordered_map<size_t, size_t> perfect_ovlp_count((*overlaps).size());
    std::unordered_map<size_t, std::vector<biosoup::Overlap>::iterator> end_filter;

    for (auto &ovlp : (*overlaps))
    {
        if (ovlp.empty())
            continue;
        std::sort(ovlp.begin(), ovlp.end(), [](biosoup::Overlap o1, biosoup::Overlap o2)
                  { return o1.score < o2.score; });
        perfect_ovlp_count[ovlp[0].lhs_id] = 0;
        for (auto &it : ovlp)
        {
            if (it.alignment == "")
            {
                continue;
            }
            // size_t begin_insertion = it.alignment.find("I");
            // size_t begin_deletion = it.alignment.find("D");
            // size_t begin_match = it.alignment.find("=");
            // size_t begin_miss = it.alignment.find("X");
            // size_t begin_indel = 0;
            // size_t min_indl = std::min(begin_insertion, begin_deletion);
            // size_t min_MM = std::min(begin_match, begin_miss);
            // if (min_indl < min_MM && min_indl != it.alignment.npos && min_MM != it.alignment.npos)
            // {
            //     begin_indel = std::stoi(it.alignment, nullptr);
            // }
            // size_t end_indel = 0;
            // if (it.alignment.at(it.alignment.size() - 1) == 'I' ||
            //     it.alignment.at(it.alignment.size() - 1) == 'D')
            // {
            //     size_t ind = it.alignment.size() - 2;
            //     while (std::isdigit(it.alignment.at(ind)))
            //     {
            //         end_indel *= 10;
            //         end_indel += it.alignment.at(ind) - '0';
            //         ind--;
            //     }
            // }
            size_t thr = (filter == filtering::DOUBLE_TRESHOLD) ? ceil((it.lhs_end - it.lhs_begin) * error_rate) : 0;
            if (it.score <= thr)
                perfect_ovlp_count[ovlp[0].lhs_id]++;
        };
    }

    switch (filter)
    {
    case filtering::NONE:
    {
        for (auto &ovlp : *overlaps)
        {
            if (ovlp.empty())
                continue;
            end_filter[ovlp[0].lhs_id] = ovlp.end();
        }
        break;
    }
    case filtering::MEDIAN:
    {
        for (auto &ovlp : (*overlaps))
        {
            if (ovlp.empty())
                continue;
            if (ovlp.size() > median_overlaps)
                end_filter[ovlp[0].lhs_id] = ovlp.begin() + median_overlaps;
            else
                end_filter[ovlp[0].lhs_id] = ovlp.end();
        }
        break;
    }
    case filtering::MEAN:
    {

        size_t round_mean = std::ceil(median_overlaps);
        for (auto &ovlp : (*overlaps))
        {
            if (ovlp.empty())
                continue;
            if (ovlp.size() > round_mean)
                end_filter[ovlp[0].lhs_id] = ovlp.begin() + round_mean;
            else
                end_filter[ovlp[0].lhs_id] = ovlp.end();
        }
        break;
    }
    case filtering::DOUBLE_PERFECT:
    {
        for (auto &ovlp : (*overlaps))
        {
            if (ovlp.empty())
                continue;
            size_t cutof = ceil(cutofparam * perfect_ovlp_count[ovlp[0].lhs_id]);
            if (ovlp.size() > cutof)
                end_filter[ovlp[0].lhs_id] = ovlp.begin() + cutof;
            else
                end_filter[ovlp[0].lhs_id] = ovlp.end();
        }
        break;
    }
    case filtering::DOUBLE_TRESHOLD:
    {
        for (auto &ovlp : (*overlaps))
        {
            if (ovlp.empty())
                continue;
            size_t cutof = std::max((size_t)(cutofparam * perfect_ovlp_count[ovlp[0].lhs_id]), median_overlaps);
            if (ovlp.size() > cutof)
                end_filter[ovlp[0].lhs_id] = ovlp.begin() + cutof;
            else
                end_filter[ovlp[0].lhs_id] = ovlp.end();
        }
    }
    }

    if (filter != filtering::NONE && print_paf)
    {
        std::cout << "printing filtered to " << out_file_filtered << std::endl;
        printPaf(overlaps.get(), end_filter, sequences, out_file_filtered);
    }

    std::cout << "setup finished overlaps " << (*overlaps).size() << std::endl;

    auto cigar_to_edlib_alignment = [](const std::string &s) -> std::string
    {
        std::string rs = "";
        std::uint64_t pos = 0;
        std::uint64_t start_pos = 0;
        std::uint64_t total_num = 0;

        if (s.empty())
        {
            return rs;
        }
        for (size_t i = 0; i < s.length(); i++)
        {
            if (std::isdigit(s[i]))
            {
                if (pos == 0)
                {
                    start_pos = i;
                }
                ++pos;
            }
            else
            {
                total_num = 0;
                for (size_t j = start_pos; j < start_pos + pos; j++)
                {
                    total_num += (s[j] - '0') * std::pow(10, (start_pos + pos) - j - 1);
                }
                pos = 0;
                switch (s[i])
                {
                case '=':
                    for (size_t j = 0; j < total_num; j++)
                    {
                        rs += '\000';
                    };
                    break;
                case 'X':
                    for (size_t j = 0; j < total_num; j++)
                    {
                        rs += '\003';
                    };
                    break;
                case 'I':
                    for (size_t j = 0; j < total_num; j++)
                    {
                        rs += '\001';
                    };
                    break;
                case 'D':
                    for (size_t j = 0; j < total_num; j++)
                    {
                        rs += '\002';
                    };
                    break;
                case 'M':
                    for (size_t j = 0; j < total_num; j++)
                    {
                        rs += '\000';
                    }
                default:
                    // rs += '\000';
                    break;
                }
            }
        }
        return rs;
    };

    auto call_snps = [&](std::uint32_t i, std::vector<biosoup::Overlap> &ovlps_final) -> void
    {
        std::uint32_t seq_inflated_len = (*sequences)[i]->inflated_len;
        std::vector<base_pile> base_pile_tmp(seq_inflated_len);
        std::vector<std::uint32_t> cov;
        double avg_mapping_quality = 0;

        // parsing alignment and creating base pile
        for (auto ovlp = ovlps_final.begin(); ovlp != end_filter[i]; ovlp++)
        {
            if (!((*ovlp).alignment.empty()))
            {
                std::uint32_t lhs_pos = (*ovlp).lhs_begin;
                std::uint32_t rhs_pos = 0;
                biosoup::NucleicAcid rhs{"",
                                         (*sequences)[(*ovlp).rhs_id]->InflateData((*ovlp).rhs_begin,
                                                                                   (*ovlp).rhs_end - (*ovlp).rhs_begin)};
                if (!(*ovlp).strand)
                    rhs.ReverseAndComplement();
                std::string rhs_tmp = rhs.InflateData();

                std::string edlib_alignment = cigar_to_edlib_alignment((*ovlp).alignment);

                for (auto &edlib_align : edlib_alignment)
                {
                    switch (edlib_align)
                    {
                    case 0:
                    case 3:
                    {
                        if (rhs_pos < rhs_tmp.size() && lhs_pos < base_pile_tmp.size())
                        {

                            switch (rhs_tmp[rhs_pos])
                            {
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
                        else // if (lhs_pos >= base_pile_tmp.size())
                        {
                            // std::cerr << "Error: Out of bounds access detected. Match" << std::endl;
                            // std::cerr << (*sequences)[ovlp.lhs_id]->name
                            //           << "\t" << (*sequences)[ovlp.lhs_id]->inflated_len // length
                            //           << "\t" << ovlp.lhs_begin
                            //           << "\t" << ovlp.lhs_end
                            //           << "\t" << (ovlp.strand ? "+" : "-")
                            //           << "\t" << (*sequences)[ovlp.rhs_id]->name
                            //           << "\t" << (*sequences)[ovlp.rhs_id]->inflated_len // length
                            //           << "\t" << ovlp.rhs_begin
                            //           << "\t" << ovlp.rhs_end
                            //           << std::endl;
                            break;
                        }
                        break;
                    }
                    case 1:
                    {
                        if (lhs_pos < base_pile_tmp.size())
                        {
                            ++base_pile_tmp[lhs_pos].i;
                        }
                        else
                        {
                            // Handle out-of-bounds access
                            // std::cerr << "Error: Out of bounds access detected. Insertion" << std::endl;
                            // std::cerr << (*sequences)[ovlp.lhs_id]->name
                            //           << "\t" << (*sequences)[ovlp.lhs_id]->inflated_len // length
                            //           << "\t" << ovlp.lhs_begin
                            //           << "\t" << ovlp.lhs_end
                            //           << "\t" << (ovlp.strand ? "+" : "-")
                            //           << "\t" << (*sequences)[ovlp.rhs_id]->name
                            //           << "\t" << (*sequences)[ovlp.rhs_id]->inflated_len // length
                            //           << "\t" << ovlp.rhs_begin
                            //           << "\t" << ovlp.rhs_end
                            //           << std::endl;
                        }
                        ++lhs_pos;
                        break;
                    }
                    case 2:
                    {
                        if (lhs_pos < base_pile_tmp.size())
                        {
                            ++base_pile_tmp[lhs_pos].d;
                        }
                        else
                        {
                            //     // Handle out-of-bounds access
                            // std::cerr << "Error: Out of bounds access detected. Deleted" << std::endl;
                            // std::cerr << (*sequences)[ovlp.lhs_id]->name
                            //           << "\t" << (*sequences)[ovlp.lhs_id]->inflated_len // length
                            //           << "\t" << ovlp.lhs_begin
                            //           << "\t" << ovlp.lhs_end
                            //           << "\t" << (ovlp.strand ? "+" : "-")
                            //           << "\t" << (*sequences)[ovlp.rhs_id]->name
                            //           << "\t" << (*sequences)[ovlp.rhs_id]->inflated_len // length
                            //           << "\t" << ovlp.rhs_begin
                            //           << "\t" << ovlp.rhs_end
                            //           << std::endl;
                        }
                        ++rhs_pos;
                        break;
                    }
                    default:
                    {
                        std::cerr << "Error: Unknown alignment detected." << std::endl;
                        break;
                        // Other cases...
                    }
                    }
                }
                avg_mapping_quality += (*ovlp).score;
            }
        }

        // calculating snps
        avg_mapping_quality /= ovlps_final.size();
        annotations[i].second = avg_mapping_quality;
        cov.reserve(base_pile_tmp.size());
        for (const auto &jt : base_pile_tmp)
        {
            if (indels)
                cov.emplace_back(jt.a + jt.c + jt.g + jt.t + jt.d + jt.i);
            else
                cov.emplace_back(jt.a + jt.c + jt.g + jt.t);
        }
        std::nth_element(cov.begin(), cov.begin() + cov.size() / 2, cov.end());
        double m = cov[cov.size() / 2] * 2. / 3.;

        std::size_t j = 0;
        for (const auto &jt : base_pile_tmp)
        {
            std::vector<double> counts;
            if (indels)
                counts = {
                    static_cast<double>(jt.a),
                    static_cast<double>(jt.c),
                    static_cast<double>(jt.g),
                    static_cast<double>(jt.t),
                    static_cast<double>(jt.d),
                    static_cast<double>(jt.i)};
            else
                counts = {
                    static_cast<double>(jt.a),
                    static_cast<double>(jt.c),
                    static_cast<double>(jt.g),
                    static_cast<double>(jt.t)};

            double sum = std::accumulate(counts.begin(), counts.end(), 0);
            if (use_frequencies)
            {
                for (auto &kt : counts)
                {
                    kt /= sum;
                };
            };

            if (sum > m)
            {
                std::size_t variants = 0;
                for (const auto &it : counts)
                {
                    if (use_frequencies)
                    {
                        if (freq_low_th < it && it < freq_high_th)
                        {
                            ++variants;
                        }
                    }
                    else
                    {
                        if (it > variant_call_th)
                        {
                            ++variants;
                        }
                    }
                }
                if (variants > 1)
                {
                    char base;
                    switch (std::distance(counts.begin(), std::max_element(counts.begin(), counts.end())))
                    {
                    case 0:
                        base = 'A';
                        break;
                    case 1:
                        base = 'C';
                        break;
                    case 2:
                        base = 'G';
                        break;
                    case 3:
                        base = 'T';
                        break;
                    case 4:
                        base = 'I';
                        break;
                    case 5:
                        base = 'D';
                        break;
                    default:
                        base = 'N';
                        break;
                    }
                    annotations[i].first.push_back(std::make_pair(static_cast<std::uint32_t>(j), base));
                }
            };
            ++j;
        };
    };

    auto threads = std::make_shared<thread_pool::ThreadPool>(64);
    std::vector<std::future<void>> futures;
    for (size_t l = 0; l < (*sequences).size(); l++)
    {
        if ((*overlaps)[l].empty())
            continue;
        // call_snps(l, (*overlaps)[l]);
        futures.emplace_back(threads->Submit(
            [&](size_t l)
            {
                call_snps(l, (*overlaps)[l]);
            },
            l));
    }
    for (auto &future : futures)
    {
        future.wait();
    }

    std::ofstream outdata;

    std::string snpOut = sample + "_" + lib + "_" + "snp_annotations.anno";
    outdata.open(snpOut);

    for (std::uint32_t i = 0; i < annotations.size(); ++i)
    {
        if (annotations[i].first.empty())
        {
            continue;
        }
        outdata << (*sequences)[i]->name;
        for (const auto &jt : annotations[i].first)
        {
            outdata << " " << jt.first << ":" << jt.second;
        }
        outdata << " mq:" << annotations[i].second;
        outdata << std::endl;
    }
    std::cout << "Calculating annotations done, writing to file" << std::endl;

    auto classify = [&](std::uint32_t i, std::vector<biosoup::Overlap> &ovlps) -> void
    {
        for (auto ovlp : ovlps)
        {
            if ((*sequences)[ovlp.rhs_id]->name == "aadf2caf-ec95-97c2-5e95-2de304acde23;chr19_PATERNAL|+strand|41613922-41719185;length=104260;error-free_length=105263;read_identity=95.748%")
            {
                int breakpoint;
                breakpoint++;
            }
            std::set<uint32_t> anno_l;

            for (const auto &anno : annotations[ovlp.lhs_id].first)
            {

                if (anno.first >= ovlp.lhs_begin && anno.first <= ovlp.lhs_end)
                {
                    anno_l.insert(anno.first);
                }
            }
            std::set<uint32_t> anno_r;
            for (const auto &anno : annotations[ovlp.rhs_id].first)
            {
                if (anno.first >= ovlp.rhs_begin && anno.first <= ovlp.rhs_end)
                {
                    anno_r.insert(anno.first);
                }
            }

            std::string edlib_alignment = cigar_to_edlib_alignment(ovlp.alignment);

            size_t lhs_it = ovlp.lhs_begin;
            size_t rhs_it = ovlp.rhs_begin;
            size_t same = 0;
            size_t different = 0;
            size_t only_on_lhs_same = 0;
            size_t only_on_rhs_same = 0;
            size_t only_on_lhs_different = 0;
            size_t only_on_rhs_different = 0;

            for (auto it : edlib_alignment)
            {
                switch (it)
                {
                case 0:
                {
                    if (anno_l.count(lhs_it) && anno_r.count(rhs_it))
                        same++;
                    if (anno_l.count(lhs_it) && !anno_r.count(rhs_it))
                        only_on_lhs_same++;
                    if (!anno_l.count(lhs_it) && anno_r.count(rhs_it))
                        only_on_rhs_same++;

                    lhs_it++;
                    rhs_it++;
                    break;
                }
                case 3:
                {
                    if (anno_l.count(lhs_it) && anno_r.count(rhs_it))
                        different++;
                    if (anno_l.count(lhs_it) && !anno_r.count(rhs_it))
                        only_on_lhs_different++;
                    if (!anno_l.count(lhs_it) && anno_r.count(rhs_it))
                        only_on_rhs_different++;

                    lhs_it++;
                    rhs_it++;
                    break;
                }
                case 1: // insertion
                {
                    lhs_it++;
                    break;
                }
                case 2: // deletion
                {
                    rhs_it++;
                    break;
                }
                }
            }

            double same_intersect = 2 * same / std::max(1UL, (only_on_lhs_same + only_on_rhs_same));
            double different_intersect = 2 * different / std::max(1UL, (only_on_lhs_different + only_on_rhs_different));

            if (!hifiasm_output)
            {
                if (same > 0.9 * (same + different) && same_intersect > 0.1)
                {
                    classified_ovlps[i].emplace_back(0, same, different, only_on_lhs_same, only_on_lhs_different, only_on_rhs_same, only_on_rhs_different);
                }
                else if (different > 0.9 * (same + different) && different_intersect > 0.1)
                {
                    classified_ovlps[i].emplace_back(1, same, different, only_on_lhs_same, only_on_lhs_different, only_on_rhs_same, only_on_rhs_different);
                }
                else if (different == 0 && same == 0 &&
                         only_on_lhs_different < no_snps_thr && only_on_lhs_same < no_snps_thr &&
                         only_on_rhs_different < no_snps_thr && only_on_rhs_same < no_snps_thr)
                {
                    classified_ovlps[i].emplace_back(3, same, different, only_on_lhs_same, only_on_lhs_different, only_on_rhs_same, only_on_rhs_different);
                }
                else
                {
                    classified_ovlps[i].emplace_back(2, same, different, only_on_lhs_same, only_on_lhs_different, only_on_rhs_same, only_on_rhs_different);
                }
            }
            else
            {
                if (different > same ||
                    (only_on_lhs_different > only_on_lhs_same &&
                     only_on_rhs_different > only_on_rhs_same))
                {
                    classified_ovlps[i].emplace_back(1, same, different, only_on_lhs_same, only_on_lhs_different, only_on_rhs_same, only_on_rhs_different);
                }
                else
                {
                    classified_ovlps[i].emplace_back(0, same, different, only_on_lhs_same, only_on_lhs_different, only_on_rhs_same, only_on_rhs_different);
                }
            }
        }
    };

    for (size_t l = 0; l < (*overlaps).size(); l++)
    {

        futures.emplace_back(threads->Submit(
            [&](size_t l)
            {
                classify(l, (*overlaps)[l]);
            },
            l));
        // classify(l, (*overlaps)[l]);
    }
    for (auto &future : futures)
    {
        future.wait();
    }

    std::cout << "Classification done, writing to files" << std::endl;

    std::ofstream out_same_haplotype;
    std::string file_sh = lib + "_" + sample + ".0.paf";
    out_same_haplotype.open(file_sh);
    std::ofstream out_different_haplotype;
    std::string file_dh = lib + "_" + sample + ".1.paf";
    out_different_haplotype.open(file_dh);
    std::ofstream out_wrong_overlaps;
    std::string wrong_overlaps = lib + "_" + sample + ".2.paf";
    out_wrong_overlaps.open(wrong_overlaps);
    std::ofstream out_no_snps;
    std::string no_snps = lib + "_" + sample + ".3.paf";
    out_no_snps.open(no_snps);

    for (size_t i = 0; i < (*overlaps).size(); i++)
    {
        for (size_t j = 0; j < (*overlaps)[i].size(); j++)
        {
            switch (std::get<0>(classified_ovlps[i][j]))
            {
            case 0:
            {
                out_same_haplotype << (*sequences)[(*overlaps)[i][j].lhs_id]->name
                                   << "\t" << (*sequences)[(*overlaps)[i][j].lhs_id]->inflated_len
                                   << "\t" << (*overlaps)[i][j].lhs_begin
                                   << "\t" << (*overlaps)[i][j].lhs_end
                                   << "\t" << ((*overlaps)[i][j].strand ? "+" : "-")
                                   << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->name
                                   << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->inflated_len
                                   << "\t" << (*overlaps)[i][j].rhs_begin
                                   << "\t" << (*overlaps)[i][j].rhs_end
                                   << "\t" << (*overlaps)[i][j].score
                                   << "\t" << 255
                                   << "\t" << std::get<1>(classified_ovlps[i][j]) << ":" << std::get<2>(classified_ovlps[i][j])
                                   << "\t" << std::get<3>(classified_ovlps[i][j]) << ":" << std::get<4>(classified_ovlps[i][j])
                                   << "\t" << std::get<5>(classified_ovlps[i][j]) << ":" << std::get<6>(classified_ovlps[i][j])
                                   << std::endl;
                break;
            }
            case 1:
            {
                out_different_haplotype << (*sequences)[(*overlaps)[i][j].lhs_id]->name
                                        << "\t" << (*sequences)[(*overlaps)[i][j].lhs_id]->inflated_len
                                        << "\t" << (*overlaps)[i][j].lhs_begin
                                        << "\t" << (*overlaps)[i][j].lhs_end
                                        << "\t" << ((*overlaps)[i][j].strand ? "+" : "-")
                                        << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->name
                                        << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->inflated_len
                                        << "\t" << (*overlaps)[i][j].rhs_begin
                                        << "\t" << (*overlaps)[i][j].rhs_end
                                        << "\t" << (*overlaps)[i][j].score
                                        << "\t" << 255
                                        << "\t" << std::get<1>(classified_ovlps[i][j]) << ":" << std::get<2>(classified_ovlps[i][j])
                                        << "\t" << std::get<3>(classified_ovlps[i][j]) << ":" << std::get<4>(classified_ovlps[i][j])
                                        << "\t" << std::get<5>(classified_ovlps[i][j]) << ":" << std::get<6>(classified_ovlps[i][j])
                                        << std::endl;
                break;
            }
            case 2:
            {
                out_wrong_overlaps << (*sequences)[(*overlaps)[i][j].lhs_id]->name
                                   << "\t" << (*sequences)[(*overlaps)[i][j].lhs_id]->inflated_len
                                   << "\t" << (*overlaps)[i][j].lhs_begin
                                   << "\t" << (*overlaps)[i][j].lhs_end
                                   << "\t" << ((*overlaps)[i][j].strand ? "+" : "-")
                                   << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->name
                                   << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->inflated_len
                                   << "\t" << (*overlaps)[i][j].rhs_begin
                                   << "\t" << (*overlaps)[i][j].rhs_end
                                   << "\t" << (*overlaps)[i][j].score
                                   << "\t" << 255
                                   << "\t" << std::get<1>(classified_ovlps[i][j]) << ":" << std::get<2>(classified_ovlps[i][j])
                                   << "\t" << std::get<3>(classified_ovlps[i][j]) << ":" << std::get<4>(classified_ovlps[i][j])
                                   << "\t" << std::get<5>(classified_ovlps[i][j]) << ":" << std::get<6>(classified_ovlps[i][j])
                                   << std::endl;
                break;
            }
            case 3:
            {
                out_no_snps << (*sequences)[(*overlaps)[i][j].lhs_id]->name
                            << "\t" << (*sequences)[(*overlaps)[i][j].lhs_id]->inflated_len
                            << "\t" << (*overlaps)[i][j].lhs_begin
                            << "\t" << (*overlaps)[i][j].lhs_end
                            << "\t" << ((*overlaps)[i][j].strand ? "+" : "-")
                            << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->name
                            << "\t" << (*sequences)[(*overlaps)[i][j].rhs_id]->inflated_len
                            << "\t" << (*overlaps)[i][j].rhs_begin
                            << "\t" << (*overlaps)[i][j].rhs_end
                            << "\t" << (*overlaps)[i][j].score
                            << "\t" << 255
                            << "\t" << std::get<1>(classified_ovlps[i][j]) << ":" << std::get<2>(classified_ovlps[i][j])
                            << "\t" << std::get<3>(classified_ovlps[i][j]) << ":" << std::get<4>(classified_ovlps[i][j])
                            << "\t" << std::get<5>(classified_ovlps[i][j]) << ":" << std::get<6>(classified_ovlps[i][j])
                            << std::endl;
            }
            }
        }
    }
    std::cout << "done" << std::endl;
    return 0;
}