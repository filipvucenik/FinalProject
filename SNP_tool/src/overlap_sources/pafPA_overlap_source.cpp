#include "overlap_source.hpp"
#include "bioparser/fasta_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

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
                sequence_id[sequence->name] = sequences.size()-1;
            }
        }

    }

    size_t find_id(std::string overlap_name){
        for(auto& it : sequences){
            if(overlap_name == it->name){
                return it->id;
            }
        }
        return -1;
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

            if(variables.size() == 12){
                tmp = variables[11];
            }else{
                tmp = "";
            }
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