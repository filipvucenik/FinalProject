#include <string>
#include <iostream>
#include <vector>

#include <filesystem>

#include "hifiasm.hpp"
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "biosoup/timer.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

std::string usage = "Usage: hifiasm <path_to_reads> <path_to_mutated_reads>\n";

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> create_parser(std::string& path){
    auto is_suffix = [] (const std::string& s, const std::string& suff) {
        return s.size() >= suff.size() && s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> ret_vector;

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

int main (int argc, char* argv[]) {
    std::filesystem::path currentDir = std::filesystem::current_path();
    std::cout << "Current working directory: " << currentDir << std::endl;
    std::cout << argc << std::endl;
    if (argc != 3) {
        std::cout << usage;
        exit(1);
    }
    std::string reads_file = argv[1];
    std::string reads_file_mutated = argv[2];

    auto target_parser = create_parser(reads_file);
    auto sequeqnce_parser = create_parser(reads_file_mutated);

    //parameters for ram
    unsigned int k = 15;
    unsigned int w = 5;
    unsigned int b = 500;
    unsigned int chain = 4;
    unsigned int matches = 100;
    unsigned int gap = 10000;
    double f = 0.001;
    auto threads = std::make_shared<thread_pool::ThreadPool>(10);

    //initializing ram minimizer_engine
    ram::MinimizerEngine minimizer_engine{
            threads,
            k,
            w,
            b,
            chain,
            matches,
            gap
    };


    while(true){
        std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
        try {
            targets = target_parser->Parse(1U << 30);
        } catch (std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            return 1;
        }

        if (targets.empty()) {
            break;
        }

        minimizer_engine.Minimize(targets.begin(), targets.end(), false);
        minimizer_engine.Filter(f);

        std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
        biosoup::NucleicAcid::num_objects = 0;

        while(true){

            std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
            try {
                sequences = sequeqnce_parser->Parse(1U << 30);
            } catch (std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
                return 1;
            }

            if (sequences.empty()) {
                break;
            }
            std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
            for (const auto& it : sequences) {
                if (it->id >= num_targets) {
                    continue;
                }
                futures.emplace_back(threads->Submit(
                        [&] (const std::unique_ptr<biosoup::NucleicAcid>& sequence)
                                -> std::vector<biosoup::Overlap> {
                            return minimizer_engine.Map(sequence, false, false, false);
                        },
                        std::ref(it)));
            }

            if (biosoup::NucleicAcid::num_objects >= num_targets) {
                break;
            }

        }
        sequeqnce_parser->Reset();
        biosoup::NucleicAcid::num_objects = num_targets;
    }

    return 0;

}