#include "overlap_source.hpp"
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ksw2.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>



class PafOverlapSource : public OverlapSource{
private:
    std::vector<std::vector<biosoup::Overlap>> overlaps;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    std::string reads;
    std::string paf_file;
    std::unordered_map<std::string , size_t> sequence_id;

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
                sequence_id[sequence->name] = sequences.size() - 1;
            }
        }

    }

    void load_paf(){
        std::string line;
        std::ifstream fileStream(paf_file);
        std::cout<<"Loading paf file"<<std::endl;
        std::string tmp;
        while(std::getline(fileStream, line)){
            std::istringstream iss(line);
            std::string v;
            std::vector<std::string> variables;
            while(std::getline(iss, v, '\t')){
                variables.push_back(v);
            }

            size_t id_l = sequence_id[variables[0]];
            size_t id_r = sequence_id[variables[5]];

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

        std::cout<<"Loaded paf file"<<std::endl;

        auto ksw2_wrapper = [&] (
                std::uint32_t i,
                const biosoup::Overlap& it,
                const std::string& lhs,
                const std::string& rhs) -> std::string {
            std::int8_t m = 3;
            std::int8_t n = -5;
            std::int8_t g = 4;
            std::int8_t e = 4;
            std::int8_t mn[25] = {
                    m, n, n, n, 0,
                    n, m, n, n, 0,
                    n, n, m, n, 0,
                    n, n, n, m, 0,
                    0, 0, 0, 0, 0
            };

            std::unordered_map<char, std::uint8_t> transform = {
                    {'A', 0}, {'a', 0},
                    {'C', 1}, {'c', 1},
                    {'G', 2}, {'g', 2},
                    {'T', 3}, {'t', 3}
            };

            auto lhs_ = new std::uint8_t[lhs.size()];
            for (std::size_t j = 0; j < lhs.size(); ++j) {
                lhs_[j] = transform[lhs[j]];
            }

            auto rhs_ = new std::uint8_t[rhs.size()];
            for (std::size_t j = 0; j < rhs.size(); ++j) {
                rhs_[j] = transform[rhs[j]];
            }

            int m_cigar = 0, n_cigar = 0;
            std::uint32_t* cigar = nullptr;

            auto score = ksw_gg2_sse(
                    nullptr,   // void *km
                    lhs.size(),  // int qlen
                    lhs_,  // const uint8_t *query
                    rhs.size(),  // int tlen
                    rhs_,  // const uint8_t *target
                    5,  // int8_t m
                    mn,  // const int8_t *mat
                    g,  // int8_t gapo
                    e,  // int8_t gape
                    500,  // int w
                    &m_cigar,  // int *m_cigar_
                    &n_cigar,  // int *n_cigar_
                    &cigar);  // uint32_t **cigar_
            std::string cigar_string = "";
            if (score > 0 && n_cigar > 0) {
                for (std::size_t j = 0; j < static_cast<std::size_t>(n_cigar); ++j) {
                    std::size_t count = cigar[j] >> 4;
                    std::size_t op = cigar[j] & 15;
                    switch (op) {
                        case 0: {  // M
                            cigar_string += std::to_string(count) + "M";
                            break;
                        }
                        case 1: {  // I
                            cigar_string += std::to_string(count) + "I";
                            break;
                        }
                        case 2: {  // D
                            cigar_string += std::to_string(count) + "D";
                            break;
                        }
                        default: break;
                    }
                }
            }
            free(cigar);
            delete[] rhs_;
            delete[] lhs_;
            return cigar_string;
        };


        auto threads = std::make_shared<thread_pool::ThreadPool>(64);
        std::vector<std::future<void>> futures;
        for(size_t i = 0; i < overlaps.size(); i++){
            futures.emplace_back(threads->Submit([&](size_t i)->void{
                for(size_t j = 0; j < overlaps[i].size(); j++){
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
                    overlaps[i][j].alignment = ksw2_wrapper(i, overlaps[i][j], lhs, rhs);
                    /*
                    if(sequences[i]->name == "read=1,reverse,position=2675766-2676093,length=327,NC_000913.3_mutated"){
                        std::cout<<overlaps[i][j].alignment<<std::endl;
                        std::cout<<"-----------------------------------------"<<std::endl;
                    }
                     */
                }
            }, i));
        }

        for(auto& future: futures){
            future.wait();
        }
        std::cout<<"pairwise alignment done"<<std::endl;

    }
public:

    std::vector<std::vector<biosoup::Overlap>>* get_overlaps() override {
        if(overlaps.empty()) {
            this->overlaps = std::vector<std::vector<biosoup::Overlap>>(sequences.size());
            load_paf();
        }

        return &overlaps;
    }
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> * get_sequences() override{
        return &sequences;
    }

    explicit PafOverlapSource(std::string& args)
    {
        std::cout<<"args: "<<args<<std::endl;
        auto spacePos = std::find(args.begin(), args.end(), ' ');
        std::cout<<"reads: "<<args.substr(0, std::distance(args.begin(), spacePos))<<std::endl;
        std::cout<<"paf_file: "<<args.substr(std::distance(args.begin(), spacePos)+1)<<std::endl;
        this->reads = args.substr(0, std::distance(args.begin(), spacePos));
        this->paf_file = args.substr(std::distance(args.begin(), spacePos)+1);
        load_sequences();
    }

};

extern "C" OverlapSource* __attribute__((visibility("default"))) create(std::string& args){
    return (OverlapSource*) new PafOverlapSource(args);
}