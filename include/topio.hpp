/** \file topio.hpp
 *
 *  \brief This file handles loading and saving topologies. Desperately need comments.
 *
 */

#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "core.hpp"
#include "global.hpp"

namespace tpp {


  /**
   *  \brief Hm. Writes topology to file, I guess?
   *
   *  \param top topology to be written
   *  \param fname output file name
   *  \param ncf ???
   *
   */
  void save_topology(const Topology & top, const char * fname, bool ncf) ;

  void save_topology_rtp(const Topology &, const char *);

  void load_topology(Topology &, const char *);

  void load_lack(Topology &, const char *);

  void check_topology(Topology &) ;

  void save_lack(const Topology &, const char *);

#if ENABLE_GAMESS_FEATURES
  void load_hessian(ublas::matrix<double>&, const char *);
#endif

  //TODO: also serialization should be included

}
#endif

