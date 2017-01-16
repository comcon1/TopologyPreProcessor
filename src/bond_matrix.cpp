#include "bond_matrix.hpp"
#include "exceptions.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

namespace tpp {
using std::string;
using std::ostringstream;
using std::endl;

using boost::format;
using boost::lexical_cast;

using OpenBabel::OBAtom;
using OpenBabel::OBMol;
using OpenBabel::OBMolBondIter;
using OpenBabel::OBMolAtomIter;
using OpenBabel::OBAtomAtomIter;
using namespace boost::numeric::ublas;

BondMatrix::BondMatrix(const OBMol &src) {
  int dim = src.NumAtoms();
  mtx = zero_matrix<bool>(dim + 1, dim + 1);
  FOR_BONDS_OF_MOL(it,const_cast<OBMol&>(src)){
  mtx(it->GetBeginAtom()->GetIdx(), it->GetEndAtom()->GetIdx()) = true;
  mtx(it->GetEndAtom()->GetIdx(), it->GetBeginAtom()->GetIdx()) = true;
}
  check();
}

void BondMatrix::print(std::ostream &os) {
  ; // echo function
}

void BondMatrix::rec_check(int _cur) {
  for (int i = 1; i < mtx.size1(); ++i) {
    if (mtx(i, _cur) && !(curset.count(i))) {
      curset.insert(i);
      rec_check(i);
    }
  }
}

bool BondMatrix::check() {
  int cur = 1;
  curset.clear();
  curset.insert(1);
  rec_check(1);
  if (curset.size() != mtx.size1() - 1) {
    ostringstream os;
    os
        << format(
            "Connected group contains %d atoms. Molecule contains %d atoms.")
            % curset.size() % (mtx.size1() - 1) << endl;
    for (std::set<int>::iterator it = curset.begin(); it != curset.end();
        ++it)
      os << (*it) << ",";
    os << " - This atoms are isolated" << endl;
    Exception e("Bad connection atoms");
    e.add("procname", "tpp::BondMatrix::check");
    e.add("error", os.str());
    throw e;
  }
  return true;
}

} // of namespace tpp
