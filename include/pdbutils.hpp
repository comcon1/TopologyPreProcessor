/**  \file pdbutils.hpp
 *
 *  \brief This file desperately needs comments.
 *
 */

#include "core.hpp"

#include <set>
#include <openbabel/mol.h>

namespace tpp {


  /**
   *
   *  \brief COMCON1 RENUMERATION ALGORITM!
   *
   *  \param verbose print additional info during execution.
   *
   */
  boost::numeric::ublas::vector<int> generate_long_tail1(OpenBabel::OBMol &, bool);

  /**
   * \brief This has something to do with mols as well, right?
   */
  AtomArray mol_renum1(OpenBabel::OBMol &, AtomArray &,
      boost::numeric::ublas::vector<int>, bool);

  namespace detail {

    void recurse_mol_renum(OpenBabel::OBMol &, AtomArray &, bub::vector<int> &,
        std::set<unsigned> &, AtomArray &, int &, int &, bool);

    void recurse_mol_scan(OpenBabel::OBMol &, std::vector<unsigned> &, unsigned,
        std::vector<unsigned> &, std::set<unsigned> &, unsigned &);

  }

}

