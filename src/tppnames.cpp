#include "tppnames.hpp"

#include <boost/lexical_cast.hpp>

using std::string;
using boost::array;
using boost::lexical_cast;

namespace tpp {


  string TTCNameGenerator::getName() {
    string a("");
    switch (instance.type) {
      case TPP_TTYPE_BON:     a = getBondName(); break;
      case TPP_TTYPE_ANG:     a = getAngleName(); break;
      case TPP_TTYPE_RBDIH:
      case TPP_TTYPE_IMPDIH:
      case TPP_TTYPE_SYMDIH:  a = getDihedralName(); break;
      //TODO: find this _unknown_ behavior
      default: a = "unknown";
    }
    return a;
  }

  TTCNameGenerator &TTCNameGenerator::set_btypes(array<string,4> vs) {
    btypes = vs;
    return (*this);
  }

  string TTCNameGenerator::getBondName() {
    string a("dfTPP_bon_");
    int first = btypes[0] > btypes[1] ? 0 : 1;
    a += btypes[first] + "_" + btypes[!first] + "_" + lexical_cast<string>(instance.dbid);
    return a;
  }

  string TTCNameGenerator::getAngleName() {
    string a("dfTPP_ang_");
    int first = btypes[0] > btypes[2] ? 0 : 2;
    int third = first ? 0 : 2;
    a += btypes[first] + "_" + btypes[1] + "_" + btypes[third];
    a += "_" + lexical_cast<string>(instance.dbid);
    return a;
  }

  string TTCNameGenerator::getDihedralName() {
    string a("dfTPP_dih_");
    if (btypes[0] > btypes[3]) 
      a += btypes[0] + "_" + btypes[1] + "_" + btypes[2] + "_" + btypes[3];
    else 
      a += btypes[3] + "_" + btypes[2] + "_" + btypes[1] + "_" + btypes[0];
    a += "_" + lexical_cast<string>(instance.dbid);
    return a;
  }

}
