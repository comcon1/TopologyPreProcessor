/** \file topio.hpp
 *
 *  \brief This file handles loading and saving topologies. Desperately need comments.
 *
 */

#ifndef TPP_TOPWRITER_H
#define TPP_TOPWRITER_H

#include "core.hpp"
#include "global.hpp"
#include "db_base.hpp"

namespace tpp {
  // predefinition
  class ForceFieldWriter;

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
        int ffID;          //!< ID of ff in database
      };

      TopologyWriter(Settings _s, const Topology &_top): settings(_s), top(_top) {;}

      /**
       *  \brief Write topology to file in GMX format.
       *
       *  \param fname output file name
       *  \param ncf no-calculate format
       *
       */
      void saveITP(const char * fname);

      /** \brief Save information about force parameters that are absent in the FF.
        *
        * \param fname output file name
        */
      void saveAbsentParametersITP(const char *fname);


      /** \brief Save pre-topology in RTP format.
        *
        * \param fname output RTP file name
        */
      void saveRTP(const char *fname);

      /** \brief Save separate force field
        *
        * \param dbsets DB connection settings
        * \param fname output ITP file name
        */
      void saveFFITP(DbBase::Settings dbsets, const char* fname);

    private:
      Settings settings;
      const Topology &top;
  };

  /** \brief Separate class writing force field as a separate file.
    *
    */
  class ForceFieldWriter: public DbBase {
    private:
        const Topology &top;
        TopologyWriter::Settings twsets;

        std::string ffDefaults;
        std::vector<std::string> atITPList; //!< list of records for [atomtypes] section

        /// Loading table of atom non-bonded parameters.
        void loadNBAtomParams();
        /// Loading FF [defaults] values
        void loadFFDefaults();
    public:

       /** \brief Constructor initializing everythin.
         *
         * \param dbopts Connection settions
         * \param twsets TopologyWriter settings
         * \param top topology object
         */
        ForceFieldWriter(DbBase::Settings dbopts_, TopologyWriter::Settings twsets_, const Topology &top_);

        /** \brief Write force field to the file.
          *
          * \param fname file name
          */
        void writeToFile(const char *fname);
  };

  //TODO: also serialization should be included

}
#endif

