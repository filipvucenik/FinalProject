#include <iostream>
#include <string>
#include "bindings/cpp/WFAligner.hpp"
#include "edlib.h"

int main()
{
    std::string q = "AAAGAAA";
    std::string t = "AAGAA";

    wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    aligner.alignEnd2End(t.c_str(), t.size(), q.c_str(), q.size());
    std::string alignment = aligner.getAlignment();

    EdlibAlignResult result = edlibAlign(
        q.c_str(), q.size(),
        t.c_str(), t.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

    std::string cigar_EDL = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
    edlibFreeAlignResult(result);

    std::cout << "EDL: " << cigar_EDL << std::endl;
    std::cout << "WFA: " << alignment << " " << aligner.getAlignmentScore() << std::endl;
    return 0;
}