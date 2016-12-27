#include "core.hpp"

namespace tpp {

  // TODO: may be remove "Input format" and "Output format"
  // TODO: may be convert to maps
  std::string in_fmt_descr(InputFormat iform) {
    switch (iform) {
    case tpp::TPP_IF_PDB:
      return "Input file format: Protein Data Bank.";
    case tpp::TPP_IF_GRO:
      return "Input file format: GROmacs structure.";
    case tpp::TPP_IF_G96:
      return "Input file format: Gromos 96 structure.";
    case tpp::TPP_IF_GAMSP:
      return "Input file format: GAMess output.";
    default:
      return "No description is available for this input format.";
    }
  }

  std::string out_fmt_descr(OutputFormat oform) {
    switch (oform) {
    case tpp::TPP_OF_PDB:
      return "Output file format: Protein Data Bank.";
    case tpp::TPP_OF_GRO:
      return "Output file format: GROmacs structure.";
    case tpp::TPP_OF_G96:
      return "Output file format: Gromos 96 structure.";
    default:
      return "No description is available for this output format.";
    }
  }

} // of namespace tpp
