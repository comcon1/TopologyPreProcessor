#ifndef TPP_ATOM_DEFINER_HEADER
#define TPP_ATOM_DEFINER_HEADER

#include "db_base.hpp"

namespace tpp {

//
//	Auxilary classses for AtomDefiner
//

// structures for 2 3 and 4 different elements
class spec2 {
private:
  int i, j;
public:
  int first() const {
    return (i < j) ? i : j;
  }
  int second() const {
    return (i < j) ? j : i;
  }
  spec2(int a, int b) :
      i(a), j(b) {
    BOOST_CHECK(a != b);
  }
};

class spec2_ {
private:
  int i, j;
public:
  int first() const {
    return i;
  }
  int second() const {
    return j;
  }
  spec2_(int a, int b) :
      i(a), j(b) {
    ;
  }
};

class spec3 {
private:
  int i, j, k;
public:
  int first() const {
    return (i < k) ? i : k;
  }
  int second() const {
    return j;
  }
  int third() const {
    return (i < k) ? k : i;
  }
  spec3(int a, int b, int c) :
      i(a), j(b), k(c) {
    BOOST_CHECK((a - b) * (a - c) * (b - c) != 0);
  }
};

class spec3_ {
private:
  int i, j, k;
public:
  int first() const {
    return i;
  }
  int second() const {
    return j;
  }
  int third() const {
    return k;
  }
  spec3_(int a, int b, int c) :
      i(a), j(b), k(c) {
    ;
  }
};

class spec4 {
private:
  int i, j, k, l;
public:
  int first() const {
    return (i < l) ? i : l;
  }
  int second() const {
    return (i < l) ? j : k;
  }
  int third() const {
    return (i < l) ? k : j;
  }
  int fourth() const {
    return (i < l) ? l : i;
  }
  spec4(int a, int b, int c, int d) :
      i(a), j(b), k(c), l(d) {
    BOOST_CHECK((a - b) * (a - c) * (a - d) * (b - c) * (b - d) * (c - d) != 0);
  }
};

class spec4_ {
private:
  int i, j, k, l;
public:
  int first() const {
    return i;
  }
  int second() const {
    return j;
  }
  int third() const {
    return k;
  }
  int fourth() const {
    return l;
  }
  spec4_(int a, int b, int c, int d) :
      i(a), j(b), k(c), l(d) {
    ;
  }
};

/**
 * \brief Class that proceed (mb produces?) atom type definition.
 */
class AtomDefiner: public DbBase {
public:
  /// Parameters for caluclations
  struct Settings{
    short ffid; /// id of current forcefield
    bool maxbonds; /// this surely means something!
    bool maxdihedrals; /// this surely means something!
    bool maxangles; /// this surely means something!
  };
  AtomDefiner(const DbBase::Settings& s1,
              const AtomDefiner::Settings& s2,
              Topology &);
  void log_scores();
  void proceed();
  void atom_align();
private:
  AtomDefiner::Settings atomSettings;
  Topology &tp;

  std::map<int, std::map<int, int> > scores; //! atom ID -> { atomtype ID -> score }

  /**
   * Function fill *scores* map according to `atom_patterns`
   * table of the database.
   */
  void smart_fit();

  //! ID -> {atomtype ID}
  std::map<int, std::set<int> > nb_suite;

  /**
   * Function matches atomtypes according only to atomic number.
   * Result is filling *nb_suite* map.
   */
  void fill_nb();

  /**
   * Summarize scores from different x_suite maps.
   * Weight coefficients of every property in atomtype definition are applied
   * here. Coef-s are defined at the top of this header in TPP_XXX defines.
   */
  void count_scores();

  /**
   * << WHAT IS THIS ?? >>
   * To spread scores found for bonded type (name-type) to all nb types
   * (uname) that corresponds to this bonded type.
   * Makes sense only if fill_bon/ang/dih are used.
   */
  void spread_atomid();

  std::map<spec2, std::set<spec2_> > bon_suite;
  std::map<spec3, std::set<spec3_> > ang_suite;
  std::map<spec4, std::set<spec4_> > dih_suite;

  void fill_bon();
  void fill_ang();
  void fill_dih();
  void print_scores(std::ostream &os);
  void smart_cgnr();

protected:
  virtual bool connect_db();

};
// of atom definer
}// of namespace tpp

#endif
