/**  \file pdbutils.hpp
 *
 *  \brief This file desperately needs comments.
 *
 */

#include "core.hpp"

#include <set>
#include <openbabel/mol.h>

using OpenBabel::OBMol;

namespace tpp {

  class Renumberer {
    protected:
      const OBMol &mol;
      bool verbosity;

      void recurseMolRenum(const AtomArray &, std::vector<unsigned> &,
          std::set<unsigned> &, AtomArray &, unsigned &, unsigned &, bool);

      /**
       * \brief Auxiliary function for recurse chain search.
       *
       */
      void recurseMolScan(std::vector<unsigned> &, unsigned,
          std::vector<unsigned> &, std::set<unsigned> &, unsigned &);

      /**
       * \brief Checks inter-consistency of the molecule and throw exceptions.
       */
      void checkMolecule();

    public:

      /** \brief Find the longest chain considering the set of parameters.
       * \param Atom set excluded from the scan
       * \param ID of starting atom. -1 for no.
       */
      std::vector<unsigned> findLongestChain(std::set<unsigned> &, int);

      /** \brief Find the longest chain in the molecule.
       *
       * Wrapper for the case of whole molecule.
       */
      std::vector<unsigned> findLongestChain();

      /**
       * \brief Renumber and rename atoms in molecule.
       * Returns array of atoms with new names and indexes
       * \param todo
       */
      AtomArray molRenumber(AtomArray &, std::vector<unsigned>, bool);

    public:

      /**
       * \brief Constructor for renumbering class
       * \param Input molecule of OpenBabel type.
       * \param Boolean verbosity level.
       */
      Renumberer(const OBMol &, bool);
  };

}

