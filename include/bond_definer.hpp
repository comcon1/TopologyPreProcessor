#ifndef TPP_BONDDEFINER_H
#define PP_BONDDEFINER_H

#include "db_base.hpp"

namespace tpp {
  /**
   * \brief A class that defines bonds between atoms, I guess?
   *
   */
  class BondDefiner: public DbBase {

    public:

      /// Calculation settings.
      struct Settings {
        int ffID;           /// forcefield ID
        bool finalize;      /// create finalized topology
        bool expanded;      /// create expanded topology
        bool verbose;       /// print additional info
      };

      BondDefiner(const DbBase::Settings& s1,
                  const BondDefiner::Settings& s2,
                  Topology &);
      ~BondDefiner();
      void bondAlign();
      void log_needed_bonds();

    private:
      BondDefiner::Settings bondSettings;
      InternalsArray bonds;
      Topology &tp;
      std::map<std::string, std::string> namemap; //! map of uname -> name
      bool genPairs;
      std::ofstream qalcfile;

      bool verbose;

      // methods
      bool connectDB() override;

      /** \brief Complete bonds exploiting DB definitions.
        */
      void fillBonds();

      /** \brief Complete angles exploting DB definitions.
        */
      void fillAngles();

      /** \brief Complete dihedrals exploiting DB definitions.
        */
      void fillDihedrals();

      /** \brief Complete extra dihedrals exploiting DB definitions [NOT IMPLEMENTED].
        */
      void fillSpecial();

      /** \brief Complete improper dihedrals exploiting DB definitions.
        */
      void fillImpropers();

      /** \brief Complete pairs according to requirements of the force field.
        */
      void fillPairs();
  };
} // of namespace tpp

#endif
