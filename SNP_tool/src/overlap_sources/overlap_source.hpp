#include "../../cmake-build-debug/_deps/biosoup-src/include/biosoup/overlap.hpp"
#include "../../cmake-build-debug/_deps/biosoup-src/include/biosoup/nucleic_acid.hpp"
#include <vector>
#include <string>
#include <memory>

class OverlapSource{
public:
    virtual  std::vector<std::vector<biosoup::Overlap>>* get_overlaps()=0;
    virtual  std::vector<std::unique_ptr<biosoup::NucleicAcid>> * get_sequences()=0;
    virtual ~OverlapSource() = default;
};

extern "C" OverlapSource* __attribute__((visibility("default"))) create(std::string& args);

