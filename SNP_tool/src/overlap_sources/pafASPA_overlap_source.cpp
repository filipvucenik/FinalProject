#include "overlap_source.hpp"
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"
#include "astarpa.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <filesystem>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

class PafOverlapSource : public OverlapSource
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
        std::ifstream fileStream(paf_file);
        std::cout << "Loading paf file" << std::endl;
        std::string tmp;

        if (!fileStream.is_open())
        {
            std::cerr << "Error opening file" << std::endl;
            std::filesystem::path currentPath = std::filesystem::current_path();
            std::cout << "Current directory: " << currentPath << std::endl;
            std::cout << "File path: " << paf_file << std::endl;
            exit(1);
        }

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
        auto astarPA_wrapper = [&](
                                   std::uint32_t i,
                                   const biosoup::Overlap &it,
                                   const std::string &lhs,
                                   const std::string &rhs) -> std::string
        {
            size_t len;
            uint8_t *cigar;
            astarpa((const uint8_t *)lhs.c_str(), lhs.size(), (const uint8_t *)rhs.c_str(), rhs.size(),
                    &cigar, &len);
            std::string cigar_string = (const char *)cigar;
            astarpa_free_cigar(cigar);

            std::string return_cigar = "";

            for (size_t i = 0; i < cigar_string.size() - 1; i++)
            {
                if (std::isalpha(cigar_string[i]) && std::isalpha(cigar_string[i + 1]))
                {
                    return_cigar += cigar_string[i] + '1';
                }
                else
                {
                    return_cigar += cigar_string[i];
                }
            }
            return_cigar += cigar_string[cigar_string.size() - 1];
            return return_cigar;
        };
        auto threads = std::make_shared<thread_pool::ThreadPool>(64);
        std::vector<std::future<void>> futures;
        for (std::size_t i = 0; i < overlaps.size(); i++)
        {
            futures.emplace_back(threads->Submit([&](std::size_t i) -> void
                                                 {
                    for(std::size_t j = 0; j < overlaps[i].size(); j++){
                        auto lhs = sequences[i]->InflateData(overlaps[i][j].lhs_begin, overlaps[i][j].lhs_end - overlaps[i][j].lhs_begin);
                        biosoup::NucleicAcid rhs_ ("", sequences[overlaps[i][j].rhs_id]->InflateData(overlaps[i][j].rhs_begin, overlaps[i][j].rhs_end - overlaps[i][j].rhs_begin));
                        if(!overlaps[i][j].strand) rhs_.ReverseAndComplement();
                        auto rhs = rhs_.InflateData();
                        /*
                        if(sequences[i]->name == "read=1,reverse,position=2675766-2676093,length=327,NC_000913.3_mutated") {
                            std::cout << lhs << std::endl;
                            std::cout << rhs << std::endl;
                        }
                         */
                        overlaps[i][j].alignment = astarPA_wrapper(i, overlaps[i][j], lhs, rhs);
                        /*
                        if(sequences[i]->name == "read=1,reverse,position=2675766-2676093,length=327,NC_000913.3_mutated"){
                            std::cout<<overlaps[i][j].alignment<<std::endl;
                            std::cout<<"-----------------------------------------"<<std::endl;
                        }
                         */
                    } }, i));
        }

        for (auto &future : futures)
        {
            future.wait();
        }
        std::cout << "pairwise alignment done" << std::endl;
    }

public:
    std::unique_ptr<std::vector<std::vector<biosoup::Overlap>>> get_overlaps() override
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

    explicit PafOverlapSource(std::string &args)
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
    return (OverlapSource *)new PafOverlapSource(args);
}

// int main()
// {
//     std::string args = "../samples/E-coli_reads_15c.fasta ../hifiasm_overlaps_15c.paf";
//     std::unique_ptr<OverlapSource> source(create(args));
//     auto overlaps = source->get_overlaps();
//     auto sequences = source->get_sequences();
//     std::cout << "done" << std::endl;
// }