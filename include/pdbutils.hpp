/**  \file pdbutils.hpp
 *
 *  \brief This file desperately needs comments.
 *
 */

#include "core.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <set>
#include <openbabel/mol.h>

namespace tpp {


/**
 * \brief This function surely does something! But what?
 */
std::pair<boost::numeric::ublas::matrix<int>,
    boost::numeric::ublas::matrix<int> > top_neig(OpenBabel::OBMol &,
    boost::numeric::ublas::matrix<int>, int, int, int, int);

/**
 *  \brief GENERATE LONGTAIL via ZOIDBERG ALGORITHM!
 *
 *  \param verbose print additional info during execution.
 *
 */
boost::numeric::ublas::vector<int> generate_long_tail(OpenBabel::OBMol &, bool verbose);

/**
 *
 *  \brief COMCON1 RENUMERATION ALGORITM!
 *
 *  \param verbose print additional info during execution.
 *
 */
boost::numeric::ublas::vector<int> generate_long_tail1(OpenBabel::OBMol &, bool verbose);

/**
 * \brief This has something to do with mols, right?
 *
 */
AtomArray mol_renum(OpenBabel::OBMol &, AtomArray &,
    boost::numeric::ublas::vector<int>,
    bool hex_flag);

/**
 * \brief This has something to do with mols as well, right?
 */
AtomArray mol_renum1(OpenBabel::OBMol &, AtomArray &,
      boost::numeric::ublas::vector<int>, bool hex_flag);
}

