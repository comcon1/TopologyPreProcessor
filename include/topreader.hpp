/** \file topreader.hpp
 *
 *  \brief This file handles loading topology data. Desperately need comments.
 *
 */

#ifndef TPP_TOPREADER_H
#define TPP_TOPREADER_H

#include "core.hpp"
#include "global.hpp"

namespace tpp {
  // @TODO: describe in details below functions

  void load_topology(Topology &, const char *);

  void load_lack(Topology &, const char *);

  void check_topology(Topology &) ;

#if ENABLE_GAMESS_FEATURES
  void load_hessian(ublas::matrix<double>&, const char *);
#endif
}

#endif
