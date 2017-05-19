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

  /** \brief Communication class for running IO procedures under needed
   * environment.
   *
   * @TODO: ... some long desc ...
   *
   */
  class TopologyIO {

    public:

      /**
       *  \brief Write topology to file in GMX format.
       *
       *  \param top topology to be written
       *  \param fname output file name
       *  \param ncf no-calculate format
       *
       */
      void saveITP(const Topology & top, const char * fname, bool ncf) ;

      /** \brief Save information about force parameters that are absent in the FF.
        *
        * \param top topology object
        * \param fname output file name
        *
        */
      void saveAbsentParametersITP(const Topology &top, const char *fname);


      /** \brief Save pre-topology in RTP format.
        *
        * \param top topology object
        * \param fname output RTP file name
        */
      void saveRTP(const Topology &top, const char *fname);
  };

  // @TODO: describe in details below functions

  void load_topology(Topology &, const char *);

  void load_lack(Topology &, const char *);

  void check_topology(Topology &) ;

#if ENABLE_GAMESS_FEATURES
  void load_hessian(ublas::matrix<double>&, const char *);
#endif

  //TODO: also serialization should be included

}
#endif

