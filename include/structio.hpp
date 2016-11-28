#ifndef TPP_STRUCTIO_H
#define TPP_STRUCTIO_H

#include "core.hpp"
#include "global.hpp"

namespace tpp {

  /** \brief Communication class for running IO procedures under needed
   * environment.
   *
   * ... some long desc ...
   *
   */
  class StructureIO {

    protected:

      /** \brief Helper function for structural recognition.
      *
      * Converts OBMol atoms to internal atoms structure.
      *
      */
      void molToAtoms(Topology &);

      bool ignoreIndexFlag;
      bool rtpoutput_file;

    public:

      /**
       * \brief The default constructor.
       *
       * \param
       * \param
       */
      StructureIO(bool ignoreIndex, bool rtpout);

      /** \brief Loading structure from std::istream
      *
      * Here is the main loading function that considers file format.
      */
      void loadFromStream(Topology &, InputFormat, std::istream *);

      /**  \brief Loading structure from file
      *
      * Just a wrapper to loadFromStream
      */
      void loadFromFile(Topology &, InputFormat, const char *);

      /** \brief Save structure to file.
      *
      * Something long..
      */
      void saveToFile(Topology &, OutputFormat, const char *);
  };

}


#endif
