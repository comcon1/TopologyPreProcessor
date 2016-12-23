#include "pdbutils.hpp"
#include "exceptions.hpp"
#include "logger.h"
#include "bond_matrix.hpp"
#include "tppnames.hpp"

#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <openbabel/mol.h>

#include <set>
#include <sstream>
#include <algorithm>

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

namespace tpp {

  void Renumberer::checkMolecule() {
    // check molecule valence of some atoms
    TPPD << "Checking bond matrix of the molecule.";

    // TODO: throw something??
    FOR_ATOMS_OF_MOL( pt, const_cast<OBMol&>(mol) ) {
      if (pt->GetValence() > 4) {
        TPPD << format("Atom %d has high valence!!") % pt->GetIdx();
        FOR_NBORS_OF_ATOM(b, &*pt ) {
          TPPD << format("--Neighbour #%d") % b->GetIdx();
        }
      } else if ( (pt->GetValence() > 1) && (pt->GetAtomicNum() == 1) ) {
        TPPD << format("Atom %d is hydrogen with high valence!!") % pt->GetIdx();
        FOR_NBORS_OF_ATOM(b, &*pt ) {
          TPPD << format("--Neighbour #%d") % b->GetIdx();
        }
      } else if ( (pt->GetValence() > 2) && (pt->GetAtomicNum() == 8) ) {
        TPPD << format("Atom %d is oxygen with high valence!!") % pt->GetIdx();
        FOR_NBORS_OF_ATOM(b, &*pt ) {
          TPPD << format("--Neighbour #%d") % b->GetIdx();
        }
      }
    }

    TPPD << "Checking bond matrix of molecule finished.";
  }

  Renumberer::Renumberer(const OBMol &_mol, bool _v): mol(_mol) {
    verbosity = _v;
    checkMolecule();
  }

  // implement function of finding the Longest chain in molecule graph
  std::vector<unsigned> Renumberer::findLongestChain(std::set<unsigned> &excluded, int startpoint) {
    std::vector<unsigned> maxtail;
    std::vector<unsigned> curtail;

    TPPD << format("Starting c1/openbabel longtail alrogithm: \n"
    "from atom no. %d, excluding %d elems. ") % startpoint % excluded.size();

    unsigned deep = 0;
    if (startpoint == -1) {
      // start recursion scan from every heavy atom
      FOR_ATOMS_OF_MOL(pA, const_cast<OBMol&>(mol) ) { // p
        if (pA->GetAtomicNum() == 1) continue;
        curtail.clear();
        curtail.push_back(pA->GetIdx());
        deep = 0;
        recurseMolScan(curtail, pA->GetIdx(), maxtail, excluded, deep);
      } // pA cycle over molecule
    } else {
      deep = 0;
      recurseMolScan(curtail, startpoint, maxtail, excluded, deep);
    }

    // output logging section
    std::ostringstream os;
    os << "Longest subchain was found:\n";
    for (auto ii: maxtail)
      os << format("%3d") % ii << std::flush;
    TPPD << os.str();

    return maxtail;
  } // end of Renumberer::findLongestChain

  // wrapper to first run of find longest chain
  std::vector<unsigned> Renumberer::findLongestChain() {
    std::set<unsigned> emptyset;
    TPPI << "Starting C1 longchain algorithm for the whole molecule.";
    auto retV = findLongestChain(emptyset, -1);
    TPPI << format("The longest chain was found: %d elements.") % retV.size();
    // logging
    return retV;
  }

  // implementation of cor renumber function
  AtomArray Renumberer::molRenumber(AtomArray &ar, std::vector<unsigned> tail, bool hexFlag) {

    TPPI << "Starting C1 renumerator alrogithm.";
  #ifdef DEBUG
    FOR_ATOMS_OF_MOL(pt, mol) {
      cout << format(" %1$3d %2$4s %3$3d\n") %pt->GetIdx() % pt->GetType() % pt->GetAtomicNum();
    }
    cout << "================================" << endl;
  #endif

    /* ATOM RENUMBERING SECTION. preparing to recursion */
    AtomArray A; // final atom array with new names and numbers
    std::set<unsigned> tailed; // set of atoms ALREADY in tails
    unsigned n = 0; // total atom counter (new index)
    unsigned h = 0; // hard atom counter

    // STARTING RECURSION
    recurseMolRenum(ar, tail, tailed, A, n, h, hexFlag);

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
      TPPD << os.str();
    } // end. output block
    return A;
  } // end mol_renum1

  // Recursion function of the renumerator algorithm
  void Renumberer::recurseMolRenum(const AtomArray &_ar, std::vector<unsigned> &_tail,
        std::set<unsigned> &_tailed, AtomArray &_A, unsigned &_n, unsigned &_h, bool hexFlag) {

      // append <tailed> array
      for (auto p: _tail)
        if (! _tailed.insert(p).second) {
          TPPE << format("Failed to insert atom no. %d") % p;
          // throwing ..
          tpp::Parameters params;
          params.add("procname", "tpp::detail::recurse_mol_scan");
          params.add("error",
            string("Error in atom numeration. See the log."));
          throw tpp::Exception("Error in recurse scan.", params);
        }
      // cycle over current tail
      for (auto p: _tail) {
        { // area of defining local variables
          OBAtom *pA = mol.GetAtom(p);
          BOOST_CHECK(_ar.count(p) == 1);
          Atom tat = *(_ar.find(p));
          _n++;
          _h++;
          tat.index = _n;
          tat.coord(0) = pA->GetX();
          tat.coord(1) = pA->GetY();
          tat.coord(2) = pA->GetZ();
          tat.ncharge = pA->GetAtomicNum();
          AtomNameGenerator ang(tat);
          tat.atom_name = ang.setNums(_h,0,hexFlag).getName();
          BOOST_CHECK(_A.insert(tat).second);
          int k = 0; // hydrogen local counter
          // arrange hydrogens
          FOR_NBORS_OF_ATOM(pQ, &*pA){
            if (pQ->GetAtomicNum() != 1) continue;
            _n++;
            k++;
            tat = * ( _ar.find(pQ->GetIdx()) );
            tat.ncharge = pQ->GetAtomicNum();
            tat.index = _n;
            tat.coord(0) = pQ->GetX();
            tat.coord(1) = pQ->GetY();
            tat.coord(2) = pQ->GetZ();
            AtomNameGenerator ang(tat);
            tat.atom_name = ang.setNums(_h,k,hexFlag).getName();
            if ( ! _A.insert(tat).second) {
              TPPE << format("Problems in renumbering hydrogen No.%d!") % pQ->GetIdx();
            }
          }
        } // end of local variables area
        // start LONGTAIL from the point and drawn deep into recursion
        while (true) {
          std::vector<unsigned> newtail = findLongestChain(_tailed, p );
          if (newtail.size() == 0) break;
          recurseMolRenum(_ar, newtail, _tailed, _A, _n, _h, hexFlag);
        }
      } // end cycle over longtail

  } // end of recurseMolRenum

  // Recursive subfunction of findLongestChain
  void Renumberer::recurseMolScan(std::vector<unsigned> &tail, unsigned cur,
        std::vector<unsigned> &maxtail, std::set<unsigned> &excluded, unsigned &deep) {
      deep++; // counter of recursion deepness
      unsigned strFwdCnt = 0; // length of non-recursive line
      if (deep < TPP_MAX_LONGTAIL_DEEP) {
          std::vector<unsigned> posNb; // possible neighbors
          unsigned curId = cur;
          while (1) {
            posNb.clear();
            FOR_NBORS_OF_ATOM(pA, &*(mol.GetAtom(curId)) ) {
              unsigned _id = pA->GetIdx();
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
                recurseMolScan(tail, jj, maxtail, excluded, deep);
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
      } // no deep limit
      else {
        TPPD << "Branching deep threshold was reached!";
        ostringstream os;
        for (auto ii: tail) os << ii << ",";
        TPPD << ( "Current: " + os.str() );
      }
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


}
