#include "pdbutils.hpp"
#include "exceptions.hpp"
#include "runtime.hpp"
#include "bond_matrix.hpp"

#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <openbabel/mol.h>

#include <set>
#include <sstream>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using std::string;
using std::ostringstream;
using std::cout;
using std::endl;

using boost::format;
using boost::lexical_cast;
using OpenBabel::OBAtom;
using OpenBabel::OBMol;
using OpenBabel::OBMolBondIter;
using OpenBabel::OBMolAtomIter;
using OpenBabel::OBAtomAtomIter;

// incapsulating some auxilary funcs
namespace {
  /**
   * Some auxiliary functions first
   * TODO: move to tppnames
   */
  string an2str(int a, string def) {
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
        return def;
    }
    return string("-");
  }

  /** Main chain private numbering function
   * TODO: move to tppnames
   */
  string mc_numerer(int id, bool hex_flag) {
    char *rrr = new char[4];
    if ((id > 255) || ( !hex_flag && (id > 99))) {
      tpp::Parameters params;
      params.add("procname", "tpp::mc_numerer");
      params.add("error",
          string("Too many atoms to number. Try to turn HEX mode. "));
      throw tpp::Exception("Bad connection atoms", params);
    }
    if (hex_flag) {
      sprintf(rrr, "%X", id);
    } else {
      sprintf(rrr, "%d", id);
    }
    return string(rrr);
  }

} // end unnamed namespace

namespace tpp {

/**
 * GENERATE LONGTAIL ALGORYTHM
 * TODO: start point!
 */
bub::vector<int> generate_long_tail1(OBMol &mol,
                                       std::set<unsigned> &excluded,
                                       int startpoint,
                                       bool verbose) {
  std::vector<unsigned> maxtail;
  std::vector<unsigned> curtail;
  int id = 0;

  runtime.log_write("Starting c1/openbabel longtail alrogithm.\n");

  unsigned deep = 0;
  if (startpoint == -1) {
    // start recursion from every heavy atom
    FOR_ATOMS_OF_MOL(pA, mol) { // p
      if (pA->GetAtomicNum() == 1) continue;
      id = pA->GetIdx();
      curtail.clear();
      curtail.push_back(id);
      deep = 0;
      detail::recurse_mol_scan(mol, curtail, id, maxtail, excluded, deep);
    } // pA cycle over molecule

  } else {
    deep = 0;
    detail::recurse_mol_scan(mol, curtail, startpoint, maxtail, excluded, deep);
  }

// generate a vector of longtail
  bub::vector<unsigned> tail_long(maxtail.size() + 1);
  for (int i = 1; i <= maxtail.size(); ++i)
    tail_long(i) = maxtail[i - 1];

  // output section
  runtime.log_write("Longest subchain is: \n");
  {
    std::ostringstream os;
    os << subrange(tail_long, 1, tail_long.size()) << std::flush;
    runtime.log_write(os.str());
    if (verbose && tail_long.size() > 1) {
      cout << "At. subchain: " << os.str() << endl;
    }
  }
  runtime.log_write("\nlongtail algroithm finished its work.\n");
  return tail_long;
}


/** Wrapper to simplier exported variant of generate_long_tail
 *
 */
bub::vector<int> generate_long_tail1(OBMol &mol, bool verbose) {
  std::set<unsigned> emptyset;
  return generate_long_tail1(mol, emptyset, -1, verbose);
}




/** Renumerate and rename atoms in molecule.
 *
 */
AtomArray mol_renum1(OBMol &mol, AtomArray &ar, bub::vector<int> tail, bool hex_flag) {

  runtime.log_write("Starting C1 renumerator alrogithm.\n");
#ifdef DEBUG
  FOR_ATOMS_OF_MOL(pt, mol) {
    cout << format(" %1$3d %2$4s %3$3d\n") %pt->GetIdx() % pt->GetType() % pt->GetAtomicNum();
  }
  cout << "================================" << endl;
#endif

  // check molecule valence of some atoms
  runtime.log_write("Checking bond matrix of molecule..\n");
  FOR_ATOMS_OF_MOL(pt, mol){
    if (pt->GetValence() > 4) {
      runtime.log_write(string("Atom ") + lexical_cast<string>(pt->GetIdx()) + " has high valence!!\n");
      FOR_NBORS_OF_ATOM(b, &*pt ) {
        runtime.log_write(string("--Neighbour #") + lexical_cast<string>(b->GetIdx()) + "\n");
      }
    } else if ( (pt->GetValence() > 1) && (pt->GetAtomicNum() == 1) ) {
      runtime.log_write(string("Atom ") + lexical_cast<string>(pt->GetIdx()) + " is hydrogen with high valence!!\n");
      FOR_NBORS_OF_ATOM(b, &*pt ) {
        runtime.log_write(string("--Neighbour #") + lexical_cast<string>(b->GetIdx()) + "\n");
      }
    } else if ( (pt->GetValence() > 2) && (pt->GetAtomicNum() == 8) ) {
      runtime.log_write(string("Atom ") + lexical_cast<string>(pt->GetIdx()) + " is oxygen with high valence!!\n");
      FOR_NBORS_OF_ATOM(b, &*pt ) {
        runtime.log_write(string("--Neighbour #") + lexical_cast<string>(b->GetIdx()) + "\n");
      }
    }
  }
  runtime.log_write("Checking bond matrix of molecule finished.\n");

  /* ATOM RENUMBERING SECTION. preparing to recursion */
  AtomArray A; // final atom array with new names and numbers
  std::set<unsigned> tailed; // set of atoms ALREADY in tails
  int n = 0; // total atom counter (new index)
  int h = 0; // hard atom counter

  // STARTING RECURSION
  detail::recurse_mol_renum(mol, ar, tail, tailed, A, n, h, hex_flag);

  { // output block
    ostringstream os;
    os << "New atom names and numbers:\n";
    for (AtomArray::iterator ait = A.begin(); ait != A.end(); ++ait) {
      os
        << format(" %1$3d %2$4s %3$3d\n") % (int) ait->index
        % ait->atom_name % (int) ait->ncharge;
    }
    os << "Table of converting numbers:\n";
    for (AtomArray::iterator ait = A.begin(); ait != A.end(); ++ait) {
      os << format("%1$3d -> %2$3d\n") % ait->oldindex % ait->index;
    }
    os << "Table of converting names:\n";
    for (AtomArray::iterator ait = A.begin(); ait != A.end(); ++ait) {
      os << format("%1$s -> %2$s\n") % ait->old_aname % ait->atom_name;
    }
    os << "COMCON1 renumeration alogrythm finished its work.\n";
    runtime.log_write(os.str());
  } // end. output block
  return A;
} // end mol_renum1


namespace detail {

  /** Recursion function of the renumerator algorithm
  *
  */
  void recurse_mol_renum(OBMol &_mol, AtomArray &_ar, bub::vector<int> &_tail,
      std::set<unsigned> &_tailed, AtomArray &_A, int &_n, int &_h, bool hex_flag) {

    // append <tailed> array
    for (int p = 1; p < _tail.size(); ++p)
      BOOST_CHECK(_tailed.insert(_tail(p)).second);
    // cycle over current tail
    for (int p = 1; p < _tail.size(); ++p) {
      { // area of defining local variables
        OBAtom *pA = _mol.GetAtom(_tail(p));
        BOOST_CHECK(_ar.count(_tail(p)) == 1);
        Atom tat = *(_ar.find(_tail(p)));
        _n++;
        _h++;
        tat.atom_name = an2str(pA->GetAtomicNum(), tat.atom_name)
          + mc_numerer(_h, hex_flag);
        tat.index = _n;
        tat.coord(0) = pA->GetX();
        tat.coord(1) = pA->GetY();
        tat.coord(2) = pA->GetZ();
        tat.ncharge = pA->GetAtomicNum();
        BOOST_CHECK(_A.insert(tat).second);
        int k = 0; // hydrogen local counter
        // arrange hydrogens
        FOR_NBORS_OF_ATOM(pQ, &*pA){
          if (pQ->GetAtomicNum() != 1) continue;
          //     BOOST_CHECK(count.insert(pQ->GetIdx()).second);
          _n++;
          k++;
          tat = * ( _ar.find(pQ->GetIdx()) );
          tat.ncharge = pQ->GetAtomicNum();
          tat.index = _n;
          tat.coord(0) = pQ->GetX();
          tat.coord(1) = pQ->GetY();
          tat.coord(2) = pQ->GetZ();
          tat.atom_name = lexical_cast<string>(k)+"H"+mc_numerer(_h, hex_flag);
          BOOST_CHECK(_A.insert(tat).second);
        }
      } // end of local variables area
      // start LONGTAIL from the point and drawn deep into recursion
      while (true) {
        bub::vector<int> newtail = generate_long_tail1(_mol, _tailed, _tail(p), hex_flag);
        if (newtail.size() == 1) break;
        recurse_mol_renum(_mol, _ar, newtail, _tailed, _A, _n, _h, hex_flag);
      }
    } // end cycle over longtail

  } // end of recurse_mol_renum

  /**
   *  Recursive subfunction of generate_long_tail
   */
  void recurse_mol_scan(OBMol &mol, std::vector<unsigned> &tail, unsigned cur,
      std::vector<unsigned> &maxtail, std::set<unsigned> &excluded, unsigned &deep) {
    deep++; // counter of recursion deepness
    if (deep > TPP_MAX_LONGTAIL_DEEP) {
      runtime.log_write("We reached the branching deep threshold!\n");
      if (tail.size() > maxtail.size()) {
        maxtail.clear();
        maxtail.insert(maxtail.begin(), tail.begin(), tail.end());
      }
      deep--;
      return;
    }
    std::vector<unsigned> posNb; // possible neighbors
    int curId = cur,
        strFwdCnt = 0; // number of non-recursive line
    while (1) {
      posNb.clear();
      FOR_NBORS_OF_ATOM(pA, &*(mol.GetAtom(curId)) ) {
        int _id = pA->GetIdx();
        if ( !( (pA->GetAtomicNum() == 1) ||
              (std::count(tail.begin(), tail.end(), _id)) ||
              (std::count(excluded.begin(), excluded.end(), _id))
              )
           )
          posNb.push_back(_id);
      } // end OB cycle
      if (posNb.size() > 1) {
        for (unsigned jj : posNb) {
          tail.push_back(jj);
          recurse_mol_scan(mol, tail, jj, maxtail, excluded, deep);
        }
        break;
      } else if (posNb.size() == 0) {
        break;
      } else {
        curId = posNb[0];
        tail.push_back(curId);
        strFwdCnt++;
        // and continue the cycle
      }
    } // cycle over strightforward sequence
    // TODO: prefer tail which contains more carbon atoms
#ifdef DEBUG
    cout << "CDB:RMS: " << std::flush;
    cout << tail;
    cout << " --]] " << std::endl;
#endif
    if (tail.size() > 0) {
        if (tail.size() > maxtail.size()) {
          maxtail.clear();
          maxtail.insert(maxtail.begin(), tail.begin(), tail.end());
        }
        if (tail.size() == strFwdCnt) // when forward chain from the begining
            tail.clear();
        else
            tail.resize(tail.size() - strFwdCnt - 1);
    }
    deep--;
  }  // end recurse_mol_scan


} // end namespace detail


}
