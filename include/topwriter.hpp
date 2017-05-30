/** \file topio.hpp
 *
 *  \brief This file handles loading and saving topologies. Desperately need comments.
 *
 */

#ifndef TPP_TOPWRITER_H
#define TPP_TOPWRITER_H

#include "core.hpp"
#include "global.hpp"

namespace tpp {

  /** \brief Communication class for running IO procedures under needed
   * environment.
   *
   * @TODO: ... some long desc ...
   *
   */
  class TopologyWriter {

    public:

      /// Parameters for caluclations
      struct Settings {
        bool ffSeparate;   //!< print separate
        bool ffFullPrint;  //!< print nonbonded & impropers also
        bool finalize;     //!< do everything to get production topology
      };

      TopologyWriter(Settings _s): settings(_s) {;}

      /**
       *  \brief Write topology to file in GMX format.
       *
       *  \param top topology to be written
       *  \param fname output file name
       *  \param ncf no-calculate format
       *
       */
      void saveITP(const Topology &top, const char * fname);

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

    private:
      Settings settings;
  };

  //TODO: also serialization should be included

}
#endif

