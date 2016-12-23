#include "tppnames.hpp"
#include "logger.h"
#include "exceptions.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <sstream>

using std::string;
using std::ostringstream;
using boost::array;
using boost::lexical_cast;
using boost::format;

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

  string AtomNameGenerator::an2str(int a) {
    switch (a) {
      case 1:
        return string("H");
      case 6:
        return string("C");
      case 7:
        return string("N");
      case 8:
        return string("O");
      case 15:
        return string("P");
      case 16:
        return string("S");
      case 17:
        return string("Cl");
      case 78:
        return string("Pt");
      default:
        TPPE << format("Unknown element no. %d") % a;
        return defaultAtomName;
    }
    return string("-");
  }

  /** \brief Main chain private numbering function
   *
   */
  string AtomNameGenerator::getName() {
    char *rrr = new char[4];
    if ( (heavyNum > 255) && !hexFlag ) {
      tpp::Parameters params;
      params.add("procname", "tpp::AtomNameGenerator::getName");
      params.add("error",
          string("Too many atoms to number. Try to turn HEX mode. "));
      throw tpp::Exception("Error in atom naming..", params);
    }
    ostringstream os;
    if (hexFlag) {
      if (instance.ncharge == 1)
        os << lightNum;
      os << format("%s%-2X") % an2str(instance.ncharge) % heavyNum;
    } else {
      if (instance.ncharge == 1)
        os << lightNum;
      os << format("%s%-2d") % an2str(instance.ncharge) % heavyNum;
    }
    return os.str();
  }

}
