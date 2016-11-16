/** \file db_scanner.hpp
 *
 *	\brief One more mysterious header that has something to do with parsing...
 *
 */

#ifndef TPP_DBSCANNER_H
#define TPP_DBSCANNER_H

#include "global.hpp"
#include "core.hpp"
#include "exceptions.hpp"

#include <mysql++/mysql++.h>
#include <set>
#include <map>


#define TPP_SMART_COEFF 10
#define TPP_ZNUC_COEFF  10000
#define TPP_BOND_COEFF  2
#define TPP_ANGLE_COEFF 4
#define TPP_DIHED_COEFF 8

namespace tpp {
  using std::set;
  // structures for 2 3 and 4 different elements
  class spec2 {
    private:
      int i,j;
    public:
      int  first() const { return (i < j) ? i : j; }
      int second() const { return (i < j) ? j : i; }
      spec2(int a, int b): i(a), j(b) { BOOST_CHECK(a != b); }
  };

  class spec2_ {
    private:
      int i,j;
    public:
      int  first() const { return i; }
      int second() const { return j; }
      spec2_(int a, int b): i(a), j(b) {;}
  };

  class spec3 {
    private:
      int i,j,k;
    public:
      int  first() const { return (i < k) ? i : k; }
      int second() const { return j; }
      int  third() const { return (i < k) ? k : i; }
      spec3(int a, int b, int c): i(a), j(b), k(c) { BOOST_CHECK((a-b)*(a-c)*(b-c) != 0); }
  };

  class spec3_ {
    private:
      int i,j,k;
    public:
      int  first() const { return i; }
      int second() const { return j; }
      int  third() const { return k; }
      spec3_(int a, int b, int c): i(a), j(b), k(c) {;}
  };

  class spec4 {
    private:
      int i,j,k,l;
    public:
      int  first() const { return (i<l) ? i : l; }
      int second() const { return (i<l) ? j : k; }
      int  third() const { return (i<l) ? k : j; }
      int fourth() const { return (i<l) ? l : i; } 
      spec4(int a, int b, int c, int d): i(a), j(b), k(c), l(d) { 
        BOOST_CHECK((a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d) != 0); 
      }
  };

  class spec4_ {
    private:
      int i,j,k,l;
    public:
      int  first() const { return i; }
      int second() const { return j; }
      int  third() const { return k; }
      int fourth() const { return l; }
      spec4_(int a, int b, int c, int d): i(a), j(b), k(c), l(d) {;}
  };

  extern bool operator <(const spec2 &, const spec2 &);
  extern bool operator <(const spec3 &, const spec3 &);
  extern bool operator <(const spec4 &, const spec4 &);

  class t_sql_exception : public Exception {
    public:
    virtual void fix_log() const; 
    t_sql_exception(const char *s, t_input_params &p): Exception(s, p) { ; }
  };

  class db_base {
    protected:
      mysqlpp::Connection *con;
      t_input_params par;
      virtual bool connect_db() throw (Exception);
    public:
      // need parameters 'host','user','dbname','password','port','ffname'
      db_base(t_input_params p) throw (Exception);
      virtual ~db_base() { delete con; }
  };
 
  /*
   * Class that accept information about DB and FF.
   */
  class db_info: public db_base {

    protected:
      virtual bool connect_db() throw (Exception);
      int ffid;
      std::string ffname;
      std::string ffdesc;
      std::string ffinclude;
      std::string ffrev;

      void getFFdata();
      void getDBdata();

    public:
      db_info(t_input_params) throw (Exception);
      int get_ffid() { return ffid; }
      std::string get_ffinclude() { return ffinclude; }
      std::string get_ffrev() { return ffrev; }
      std::string get_statistics();

  };

  /*
   * Class that proceed atom type definition.
   */
  class atom_definer: public db_base {
   private:
    Topology &tp;

    /*
     * atom ID -> 
     *          { atomtype ID -> score }
     */
    std::map<int, std::map<int, int> > scores;

    /*
     * Function fill *scores* map according to `atom_patterns`
     * table of the database.
     */
    void smart_fit() throw (Exception);
    
    /*
     * ID -> {atomtype ID}
     */ 
    std::map<int, std::set<int> > nb_suite;

    /*
     * Function matches atomtypes according only to atomic number.
     * Result is filling *nb_suite* map.
     */
    void fill_nb() throw (Exception);
    
    /*
     * Summarize scores from different x_suite maps.
     * Weight coefficients of every property in atomtype definition are applied
     * here. Coef-s are defined at the top of this header in TPP_XXX defines.
     */
    void count_scores() throw (Exception);

    /*
     * << WHAT IS THIS ?? >>
     * To spread scores found for bonded type (name-type) to all nb types
     * (uname) that corresponds to this bonded type. 
     * Makes sense only if fill_bon/ang/dih are used.
     */
    void spread_atomid() throw (Exception);

    std::map<spec2, std::set<spec2_> > bon_suite;
    std::map<spec3, std::set<spec3_> > ang_suite;
    std::map<spec4, std::set<spec4_> > dih_suite;

    short  ffid; // id of current forcefield
    void fill_bon() throw (Exception);
    void fill_ang() throw (Exception);
    void fill_dih() throw (Exception);
    void print_scores(std::ostream &os);
    void smart_cgnr() throw (Exception);

   protected:
    virtual bool connect_db() throw (Exception);

   public:

    atom_definer(t_input_params, Topology &) throw (Exception);
    void log_scores();
    void proceed() throw (Exception);
    void atom_align() throw (Exception);
  };

  class bond_definer: public db_base {
    private:
      InternalsArray bonds;
      Topology &tp;
      std::map<std::string, std::string> namemap; // map of uname -> name ))
      short ffid;
      bool  genpairs;
      std::ofstream qalcfile;

      // methods
      virtual bool connect_db() throw (Exception);
      void fill_bonds() throw (Exception);
      void fill_angles() throw (Exception);
      void fill_dihedrals() throw (Exception);
      void fill_special() throw (Exception);
      void fill_impropers() throw (Exception, DbException);
      void fill_pairs() throw (Exception);
    public:
      bond_definer(t_input_params, Topology &) throw (Exception);
      virtual ~bond_definer(); 
      void bond_align() throw (Exception);
      void log_needed_bonds();
  };

}

#endif

