#ifndef OVERLAP_SOURCE_HPP
#define OVERLAP_SOURCE_HPP

#include <vector>
#include <string>
#include <memory>

#include "biosoup/overlap.hpp"
#include "biosoup/nucleic_acid.hpp"

class OverlapSource
{
public:
    virtual std::vector<std::vector<biosoup::Overlap>> *get_overlaps() = 0;
    virtual std::vector<std::unique_ptr<biosoup::NucleicAcid>> *get_sequences() = 0;
    virtual ~OverlapSource() = default;
};

// extern "C" OverlapSource *__attribute__((visibility("default"))) create(std::string &args);
#endif // OVERLAP_SOURCE_HPP