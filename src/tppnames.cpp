#include "tppnames.hpp"
#include "logger.hpp"
#include "exceptions.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <sstream>
#include <cctype>

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

  /**
   *  \brief Conversion to base36 digit
   *
   *  Function is local and not used outside the module.
   */
  std::string to_base36(unsigned int val) {
      static std::string base36 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
      std::string result;
      result.reserve(14);
      do {
          result = base36[val % 36] + result;
      } while (val /= 36);
      return result;
  }

  /**
   *  \brief Namer initialization
   *
   */
  AtomNameGenerator &AtomNameGenerator::setNums(unsigned hv_, unsigned lt_, bool flag_) {

    heavyNum = hv_;
    lightNum = lt_;
    b36Flag = flag_;

    if ( (heavyNum > 99) && !b36Flag ) {
      tpp::Exception e("Error in atom naming.");
      e.add("procname", "tpp::AtomNameGenerator::getName");
      e.add("error", "Too many atoms to number. Try to turn Base36 mode. ");
      throw e;
    }

    if ( heavyNum > 1295 ) {
      tpp::Exception e("Error in atom naming.");
      e.add("procname", "tpp::AtomNameGenerator::getName");
      e.add("error", "Too many atoms for unique naming. Try to decrease the residue. ");
      throw e;
    }

    return *this;
  }


  /** \brief Main chain private numbering function
   *
   */
  string AtomNameGenerator::getName() {
    char *rrr = new char[4];

    ostringstream os;
    if (b36Flag) {
      if (instance.ncharge == 1)
        os << format("%1d%1s%s") % lightNum % an2str(instance.ncharge) % to_base36(heavyNum);
      else
        os << format("%2s%s")   % an2str(instance.ncharge) % to_base36(heavyNum);
    } else {
      if (instance.ncharge == 1)
        os << format("%1d%1s%-2d") % lightNum % an2str(instance.ncharge) % heavyNum;
      else
        os << format("%2s%-2d") % an2str(instance.ncharge) % heavyNum;
    }
    return os.str();
  } // end ANG::getName

  /** \brief Guess residue name from filename
    */
  string ResidueNameGenerator::getName() {
    string rsn;
    for (auto i: instance) {
      if ( std::isalnum(i) ) {
        rsn += toupper(i);
      }
    }
    if (rsn.length() < 3) {
      if (rsn.length() == 2) rsn += "X";
      else if (rsn.length() == 1) rsn += "XX";
      else { rsn = "XXX"; }
    } else {
      rsn = rsn.substr(0,3);
    }
    return rsn;
  }

}
