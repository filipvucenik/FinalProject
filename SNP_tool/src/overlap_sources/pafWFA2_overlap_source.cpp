#include "overlap_source.hpp"
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"
#include "bindings/cpp/WFAligner.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <filesystem>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

enum class operation
{
    MATCH,
    MISMATCH,
    INSERTION,
    DELETION,
    UNDEFINED
};

class PafWFA2OverlapSource : public OverlapSource
{
private:
    std::vector<std::vector<biosoup::Overlap>> overlaps;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    std::string reads;
    std::string paf_file;
    std::unordered_map<std::string, std::size_t> sequence_id;

    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> create_parser(std::string &path)
    {
        auto is_suffix = [](const std::string &s, const std::string &suff)
        {
            return s.size() >= suff.size() && s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
        };

        if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
            is_suffix(path, ".fna") || is_suffix(path, ".fna.gz") ||
            is_suffix(path, ".fa") || is_suffix(path, ".fa.gz"))
        {
            try
            {
                return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);
            }
            catch (const std::invalid_argument &exception)
            {
                std::cerr << exception.what() << std::endl;
            }
        }
        return nullptr;
    }
    std::string op_and_num_to_str(operation op, std::size_t num, size_t *edit_distance)
    {
        switch (op)
        {
        case operation::MATCH:
            return std::to_string(num) + "=";
        case operation::MISMATCH:
            (*edit_distance) += num;
            return std::to_string(num) + "X";
        case operation::DELETION:
            (*edit_distance) += num;
            return std::to_string(num) + "D";
        case operation::INSERTION:
            (*edit_distance) += num;
            return std::to_string(num) + "I";
        default:
            return "";
        }
    }

    void load_sequences()
    {
        auto parser = create_parser(reads);
        while (true)
        {
            std::vector<std::unique_ptr<biosoup::NucleicAcid>> buffer;
            try
            {
                buffer = parser->Parse(1U << 30);
            }
            catch (std::invalid_argument &exception)
            {
                std::cerr << exception.what() << std::endl;
                exit(1);
            }

            if (buffer.empty())
            {
                break;
            }
            sequences.reserve(sequences.size() + buffer.size());
            for (const auto &sequence : buffer)
            {
                sequences.push_back(std::make_unique<biosoup::NucleicAcid>(*sequence));
                sequence_id[sequence->name] = sequences.size() - 1;
            }
        }
    }

    std::size_t find_id(std::string overlap_name)
    {
        for (auto &it : sequences)
        {
            if (overlap_name == it->name)
            {
                return it->id;
            }
        }
        return -1;
    }

    void load_paf()
    {
        std::string line;
        if (!std::filesystem::exists(paf_file))
        {
            std::cerr << "File does not exist: " << paf_file << std::endl;
        }
        std::ifstream fileStream(paf_file);
        if (!fileStream.is_open())
        {
            std::cerr << "Error opening file" << std::endl;
            std::filesystem::path currentPath = std::filesystem::current_path();
            std::cout << "Current directory: " << currentPath << std::endl;
            std::cout << "File path: " << paf_file << std::endl;
            exit(1);
        }

        std::cout << "Loading paf file" << std::endl;
        std::string tmp;
        while (std::getline(fileStream, line))
        {
            std::istringstream iss(line);
            std::string v;
            std::vector<std::string> variables;
            while (std::getline(iss, v, '\t'))
            {
                variables.push_back(v);
            }

            std::size_t id_l = sequence_id[variables[0]];
            std::size_t id_r = sequence_id[variables[5]];

            overlaps[id_l].emplace_back(id_l,
                                        std::stoi(variables[2]),
                                        std::stoi(variables[3]),
                                        id_r,
                                        std::stoi(variables[7]),
                                        std::stoi(variables[8]),
                                        255,
                                        tmp,
                                        variables[4] == "+");
        }

        std::cout << "Loaded paf file" << std::endl;
        auto wfa2_wrapper = [&](
                                std::uint32_t i,
                                const biosoup::Overlap &it,
                                const std::string &lhs,
                                const std::string &rhs) -> std::pair<std::string, int>
        {
            wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
            aligner.alignEnd2End(rhs.c_str(), rhs.size(), lhs.c_str(), lhs.size());
            std::string alignment = aligner.getAlignment();
            size_t edit_distance = 0;
            std::size_t lhs_i = 0;
            std::size_t rhs_i = 0;
            std::string cigar_string = "";
            operation current_op = operation::UNDEFINED;
            std::size_t inarw_count = 0;
            for (std::size_t j = 0; j < alignment.size(); j++)
            {
                if (alignment[j] == 'M')
                {
                    if (current_op != operation::MATCH)
                    {
                        cigar_string += op_and_num_to_str(current_op, inarw_count, &edit_distance);
                        current_op = operation::MATCH;
                        inarw_count = 0;
                    }
                    lhs_i++;
                    rhs_i++;
                }
                else if (alignment[j] == 'X')
                {
                    if (current_op != operation::MISMATCH)
                    {
                        cigar_string += op_and_num_to_str(current_op, inarw_count, &edit_distance);
                        current_op = operation::MISMATCH;
                        inarw_count = 0;
                    }
                    lhs_i++;
                    rhs_i++;
                }
                else if (alignment[j] == 'I')
                {
                    if (current_op != operation::INSERTION)
                    {
                        cigar_string += op_and_num_to_str(current_op, inarw_count, &edit_distance);
                        current_op = operation::INSERTION;
                        inarw_count = 0;
                    }
                    lhs_i++;
                }
                else if (alignment[j] == 'D')
                {
                    if (current_op != operation::DELETION)
                    {
                        cigar_string += op_and_num_to_str(current_op, inarw_count, &edit_distance);
                        current_op = operation::DELETION;
                        inarw_count = 0;
                    }
                    rhs_i++;
                }
                inarw_count++;
            }
            cigar_string += op_and_num_to_str(current_op, inarw_count, &edit_distance);
            return std::make_pair(cigar_string, edit_distance);
        };

        auto threads = std::make_shared<thread_pool::ThreadPool>(32);
        std::vector<std::future<void>> futures;

        wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);

        for (std::size_t i = 0; i < overlaps.size(); i++)
        {
            futures.emplace_back(threads->Submit([&](std::size_t i) -> void
                                                 {
                    for(std::size_t j = 0; j < overlaps[i].size(); j++){

                        auto lhs = sequences[i]->InflateData(overlaps[i][j].lhs_begin, overlaps[i][j].lhs_end - overlaps[i][j].lhs_begin);
                        biosoup::NucleicAcid rhs_ ("", sequences[overlaps[i][j].rhs_id]->InflateData(overlaps[i][j].rhs_begin, overlaps[i][j].rhs_end - overlaps[i][j].rhs_begin));
                        if(!overlaps[i][j].strand) rhs_.ReverseAndComplement();
                        auto rhs = rhs_.InflateData();

                        std::pair<std::string, int> res = wfa2_wrapper(i, overlaps[i][j], lhs, rhs);
                        overlaps[i][j].alignment = res.first;
                        overlaps[i][j].score = res.second;
                        // if(overlaps[i][j].lhs_id == 2 && overlaps[i][j].rhs_id == 77){
                        //     std::cout << "lhs: " << sequences[i]->name << std::endl << lhs << std::endl;
                        //     std::cout << "rhs: " << sequences[overlaps[i][j].rhs_id]->name << std::endl <<rhs<< std::endl;
                        //     std::cout << "alignment: " << overlaps[i][j].alignment << std::endl;
                        // }

                    } }, i));

            // single thread
            // for (std::size_t j = 0; j < overlaps[i].size(); j++)
            // {
            //     auto lhs = sequences[i]->InflateData(overlaps[i][j].lhs_begin, overlaps[i][j].lhs_end - overlaps[i][j].lhs_begin);
            //     biosoup::NucleicAcid rhs_("", sequences[overlaps[i][j].rhs_id]->InflateData(overlaps[i][j].rhs_begin, overlaps[i][j].rhs_end - overlaps[i][j].rhs_begin));
            //     if (!overlaps[i][j].strand)
            //         rhs_.ReverseAndComplement();
            //     auto rhs = rhs_.InflateData();

            //     overlaps[i][j].alignment = wfa2_wrapperSGT(i, overlaps[i][j], lhs, rhs, aligner);
            // }
        }

        for (auto &future : futures)
        {
            future.wait();
        }
        std::cout << "pairwise alignment done" << std::endl;
    }

public:
    std::unique_ptr<std::vector<std::vector<biosoup::Overlap>>>
    get_overlaps() override
    {
        if (overlaps.empty())
        {
            this->overlaps = std::vector<std::vector<biosoup::Overlap>>(sequences.size());
            load_paf();
        }

        return std::make_unique<std::vector<std::vector<biosoup::Overlap>>>(std::move(overlaps));
    }
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> *get_sequences() override
    {
        return &sequences;
    }

    explicit PafWFA2OverlapSource(std::string &args)
    {
        std::cout << "args: " << args << std::endl;
        auto spacePos = std::find(args.begin(), args.end(), ' ');
        std::cout << "reads: " << args.substr(0, std::distance(args.begin(), spacePos)) << std::endl;
        std::cout << "paf_file: " << args.substr(std::distance(args.begin(), spacePos) + 1) << std::endl;
        this->reads = args.substr(0, std::distance(args.begin(), spacePos));
        this->paf_file = args.substr(std::distance(args.begin(), spacePos) + 1);
        load_sequences();
    }
};

extern "C" OverlapSource *__attribute__((visibility("default"))) create(std::string &args)
{
    return (OverlapSource *)new PafWFA2OverlapSource(args);
}
