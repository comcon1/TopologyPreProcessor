#ifndef TPP_ATOM_DEFINER_HEADER
#define TPP_ATOM_DEFINER_HEADER

#include "db_base.hpp"

#include <stdexcept>

namespace tpp {
  using namespace boost::multi_index;

  /**
    *	\bief Auxilary classes for AtomDefiner
    *
    * Structures for 2 3 and 4 different elements and comparison operators
    * under these objects.
    */
  namespace detail {

    /** \brief Pair with first less than second
      */
    template<typename ElemType>
    class Spec2 {
      private:
        ElemType i, j;
      public:
        ElemType first() const  { return (i < j) ? i : j; }
        ElemType second() const { return (i < j) ? j : i; }
        Spec2(ElemType a, ElemType b): i(a), j(b) { 
            if (a == b) throw std::logic_error("violation spec2 type"); }
    };

    /** \brief Simple pair<int,int>
      */
    template<typename ElemType>
    class Spec2_ {
      private:
        ElemType i, j;
      public:
        ElemType first() const { return i; }
        ElemType second() const { return j; }
        Spec2_(ElemType a, ElemType b): i(a), j(b) { ; }
    };

    /** \brief Trinity of numbers a,b,c with orientation: a < c != b
      *
      */
    template<typename ElemType>
    class Spec3 {
      private:
        ElemType i, j, k;
      public:
        ElemType first() const { return (i < k) ? i : k; }
        ElemType second() const { return j; }
        ElemType third() const { return (i < k) ? k : i; }
        Spec3(ElemType a, ElemType b, ElemType c): i(a), j(b), k(c) {
          if (! ( (a != b) && (a != c) && (b != c) ) )
              throw std::logic_error("violation spec3 type");
        }
    };

    /** \brief Simple trinity of numbers.
      */
    template<typename ElemType>
    class Spec3_ {
      private:
        ElemType i, j, k;
      public:
        ElemType first() const { return i; }
        ElemType second() const { return j; }
        ElemType third() const { return k; }
        Spec3_(ElemType a, ElemType b, ElemType c): i(a), j(b), k(c) { ; }
    };

    /** \brief Four numbers a,b,c,d with orientation a > d != b != c
      */
    template<typename ElemType>
    class Spec4 {
      private:
        ElemType i, j, k, l;
      public:
        ElemType first()  const { return (i < l) ? i : l; }
        ElemType second() const { return (i < l) ? j : k; }
        ElemType third()  const { return (i < l) ? k : j; }
        ElemType fourth() const { return (i < l) ? l : i; }
        Spec4(ElemType a, ElemType b, ElemType c, ElemType d) :
            i(a), j(b), k(c), l(d) {
          if ( !( (a != b) && (a != c) && (a != d) && (b != c) && (b != d) && (c != d) ))
              throw std::logic_error("violation spec4 type");
        }
    };

    /** \brief Simple quad of numbers.
      */
    template<typename ElemType>
    class Spec4_ {
    private:
      ElemType i, j, k, l;
    public:
      ElemType first()  const { return i; }
      ElemType second() const { return j; }
      ElemType third()  const { return k; }
      ElemType fourth() const { return l; }
      Spec4_(ElemType a, ElemType b, ElemType c, ElemType d):
          i(a), j(b), k(c), l(d) { ;  }
    };

    template<typename ElemType>
    bool operator <(const Spec2<ElemType> &a, const Spec2<ElemType> &b) {
      return (a.first() < b.first())
          || ((a.first() == b.first()) && (a.second() < b.second()));
    }

    template<typename ElemType>
    bool operator <(const Spec2_<ElemType> &a, const Spec2_<ElemType> &b) {
      return (a.first() < b.first())
          || ((a.first() == b.first()) && (a.second() < b.second()));
    }

    template<typename ElemType>
    bool operator <(const Spec3<ElemType> &a, const Spec3<ElemType> &b) {
      return (a.first() < b.first())
          || ((a.first() == b.first())
              && (Spec2<ElemType>(a.second(), a.third())
                  < Spec2<ElemType>(b.second(), b.third())));
    }

    template<typename ElemType>
    bool operator <(const Spec3_<ElemType> &a, const Spec3_<ElemType> &b) {
      return (a.first() < b.first())
          || ((a.first() == b.first())
              && (Spec2_<ElemType>(a.second(), a.third())
                  < Spec2_<ElemType>(b.second(), b.third())));
    }

    template<typename ElemType>
    bool operator <(const Spec4<ElemType> &a, const Spec4<ElemType> &b) {
      return (a.first() < b.first())
          || ((a.first() == b.first())
              && (Spec3<ElemType>(a.second(), a.third(), a.fourth())
                  < Spec3<ElemType>(b.second(), b.third(), b.fourth())));
    }

    template<typename ElemType>
    bool operator <(const Spec4_<ElemType> &a, const Spec4_<ElemType> &b) {
      return (a.first() < b.first())
          || ((a.first() == b.first())
              && (Spec3_<ElemType>(a.second(), a.third(), a.fourth())
                  < Spec3_<ElemType>(b.second(), b.third(), b.fourth())));
    }

  } // end namespace detail

  /**
   * \brief Class that performs atom type definition.
   */
  class AtomDefiner: public DbBase {
    public:
      /// Parameters for caluclations
      struct Settings {
        short ffID;        //!< id of current forcefield
        bool verbose;      //!< if the output is verbosed
        bool maxbonds;     //!< try to maximize count of known bonds [EXPERIMENTAL]
        bool maxdihedrals; //!< try to maximize count of known dihedrals [EXPERIMENTAL]
        bool maxangles;    //!< try to maximize count of known angles [EXPERIMENTAL]
      };

      /// Structure for atomtype table loading
      struct tempstruct_t {
        int id;
        std::string type;
        std::string type2;
        double charge;
        double mass;
        std::string comment;
      }; // @TODO: rename properly

      /// Type for atomtype table loading
      typedef multi_index_container<tempstruct_t,
        indexed_by<
            ordered_unique<member<tempstruct_t, int,
              &tempstruct_t::id> > > > AtomMapper; // @TODO: rename properly

      /** \brief Constructor defines initial parameters of atom definition
        */
      AtomDefiner(const DbBase::Settings&,
                  const AtomDefiner::Settings&,
                  Topology &);

      /** \brief Print current scores into log-file
        */
      void logScores();

      /** \brief The first stage procedure: filling scores tables.
        */
      void proceed();

      /** \brief The second stage procedure: atomtype attribution.
        */
      void atomAlign();

    private:
      AtomDefiner::Settings atomSettings;
      Topology &tp;

      //! atom ID -> { atomtype ID -> score }
      std::map<int, std::map<int, int> > scores;

      /** Generate initial *scores* map, filling it with zero-sets
        */
      void scoresZeroFill();

      /**
       * Function fill *scores* map according to `atom_patterns`
       * table of the database.
       */
      void smartFit();

      /**
       * \brief Function matches atomtypes according only to atomic number.
       * Result is filling *nbSuite* map.
       */
      void fillNB();

      std::map<int, std::set<std::string> > nbSuite;      //!< atom ID -> {atom valence type}
      std::map<int, std::set<std::string> > mapAItoAVT;   //! atomuc.num -> {A.V.T. set}

      std::string getSQLSet(int);                      //! get a.v.t. set in SQL form by atomic num

      std::map<detail::Spec2<int>, std::set<detail::Spec2_<std::string> > > bondSuite;   //!< (atom ID, atom ID) -> { (a.v.t, a.v.t.) }
      std::map<detail::Spec3<int>, std::set<detail::Spec3_<std::string> > > angleSuite;  //!< (atom ID, atom ID, atom ID) -> { (a.v.t, a.v.t., a.v.t.) }
      std::map<detail::Spec4<int>, std::set<detail::Spec4_<std::string> > > dihdSuite;   //!< (atom ID, atom ID, atom ID, atom ID) -> { (a.v.t, a.v.t., a.v.t., a.v.t.) }

      std::map<int, std::map<std::string, int> > avtScores;                 //!< SCORE for valence type (name) of an atom

      /**
       * \brief Function matches bondtypes according to existing bonds in DB.
       * Result is filling *bondSuite* map.
       */
      void fillBonds();

      ///TODO: REFACTOR AS fillBonds
      void fillAngles();
      void fillDihs();

      /**
       * Summarize scores from different x_suite maps.
       * Weight coefficients of every property in atomtype definition are applied
       * here. Coef-s are defined at the top of this header in TPP_XXX defines.
       */
      void countAVTScores();

      /**
       * \brief Apply scores calculated in fill_bon/ang/dih to global atomtype scores.
       *
       * To spread scores found for bonded type (name-type) to all nb types
       * (uname) that corresponds to this bonded type.
       * Makes sense only if fill_bon/ang/dih are used.
       */
      void convertAVTtoScores();


      /** \brief Print scores into stream
        */
      void printScores(std::ostream &os);

      /** \brief Charge group definition
        */
      void smartCgnr();

    protected:
      AtomMapper atom_mapper; //!< map of atomtypes loaded completely from DB

      /** \brief Now this function is not overriden.
        * See parent function for details.
        */
      virtual bool connectDB();

      /** \brief Prints useful information about SMART fit procedure results.
        *
        *  According to this output user can decide if SMART fit was performed
        * good or not. All fit results are summarized in a single table and fill_*
        * scores (AVT scores) are not considered here.
        */
      void printSmartFitStats(std::map<int, std::map<int,int> > &_sfscores,
        std::map<int, std::map<int, std::string> > &_sfsmarts);

      /** \brief Checks if atoms in tp.atoms and atoms in tp.mol are consistent.
        */
      void checkAtomlistConsistency();
  };
  // of atom definer
}// of namespace tpp

#endif
