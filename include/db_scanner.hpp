#ifndef TPP_DBSCANNER_H
#define TPP_DBSCANNER_H

#include "global.hpp"
#include "mysql++.h"
#include <set>


#define TPP_SMART_COEFF 100
#define  TPP_ZNUC_COEFF 10000
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
      spec4(int a, int b, int c, int d): i(a), j(b), k(c), l(d) { BOOST_CHECK((a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d) != 0); }
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

  class t_sql_exception : public t_exception {
    public:
    virtual void fix_log() const {
      std::ostringstream os; 
      os << "TPP catched exception!\n";
      os << format("***** from %1% -> %2%\n") % PARAM_READ(pars, "classname") % PARAM_READ(pars, "procname");
      os << "***** " << mesg << endl;
      os << "***** MYSQL: " << PARAM_READ(pars, "sql_error") << endl;
      runtime.log_write(os.str());
    }
    t_sql_exception(const char *s, t_input_params &p) { mesg = s; pars = p; }
  };
 
  class atom_definer {
   private:
    // super)
    t_input_params par;
    t_topology &tp;
    // mysql variables
    mysqlpp::Connection *con;

    map<int, map<int, int> > scores;
    map<int, set<int> > nb_suite;       // suite on znuc
    map<spec2, set<spec2_> > bon_suite;
    map<spec3, set<spec3_> > ang_suite;
    map<spec4, set<spec4_> > dih_suite;

    short  ffid; // id of current forcefield
    bool connect_db() throw (t_exception);
    void fill_nb() throw (t_exception);
    void fill_bon() throw (t_exception);
    void fill_ang() throw (t_exception);
    void fill_dih() throw (t_exception);
    void spread_atomid() throw (t_exception);
    void smart_fit() throw (t_exception);
    void count_scores() throw (t_exception);
    void print_scores(std::ostream &os);

   public:

    atom_definer(t_input_params, t_topology &) throw (t_exception);
    // need parameters 'host','user','dbname','password','port','ffname'
    ~atom_definer() 
      { delete con; }
    void log_scores();
    void proceed() throw (t_exception);
    void atom_align() throw (t_exception);
  };

  class bond_definer {
    private:
      t_internals_array bonds;
      mysqlpp::Connection *con;
      t_topology &tp;
      t_input_params par;
      map<string, string> namemap; // map of uname -> name ))
      short ffid;
      std::ofstream qalcfile;

      // methods
      bool connect_db() throw (t_exception);
      void fill_bonds() throw (t_exception);
      void fill_angles() throw (t_exception);
      void fill_dihedrals() throw (t_exception);
      void fill_special() throw (t_exception);
    public:
      bond_definer(t_input_params, t_topology &) throw (t_exception);
      ~bond_definer() {        
        delete con;
        qalcfile.close();
      }
      void bond_align() throw (t_exception);
      void log_needed_bonds();
  };
}

#endif

