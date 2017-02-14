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
        std::string ffName; /// forcefield name
        bool noqalculate;   /// What an interesting word, 'noqalculate'
        bool verbose;       /// print additional info
      };

      BondDefiner(const DbBase::Settings& s1,
                  const BondDefiner::Settings& s2,
                  Topology &);
      ~BondDefiner();
      void bond_align();
      void log_needed_bonds();

    private:
      BondDefiner::Settings bondSettings;
      InternalsArray bonds;
      Topology &tp;
      std::map<std::string, std::string> namemap; //! map of uname -> name
      short ffID;    // why not in settings?
      bool genPairs;
      std::ofstream qalcfile;

      bool verbose;

      // methods
      bool connectDB() override;

      /** \brief Complete bonds exploiting DB definitions.
        */
      void fillBonds();

      void fill_angles();
      void fill_dihedrals();
      void fill_special();
      void fill_impropers();
      void fill_pairs();
  };
} // of namespace tpp

#endif
