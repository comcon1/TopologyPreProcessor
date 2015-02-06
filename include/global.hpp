#ifndef TPP_GLOBAL_H
#define TPP_GLOBAL_H

#define PROG_NAME "TPP-0.6.1 (testing version)"
#define VERSION "0.6.1"

// system
#include <unistd.h>
// STL
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iomanip>
#include <map>
#include <memory>
// smart IO features
#include <boost/lexical_cast.hpp>
#include <boost/cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
// smart serialization features
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
// linear algebra
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
// smart containers
#include <boost/array.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
// Open Babel
#include <openbabel/mol.h>
// internal class for comparison of float/double/vector
#include <floatcomp.h> 


using std::string;
using std::pair;
using std::map;
using std::vector;
using std::string;
using std::fstream;
using std::ostream;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;
using std::auto_ptr;
using std::ios;
using std::ostringstream;

#define TPP_GMX_EXTERNAL
#define TPP_MAX_ATOM_NUM 200
#define TPP_INDEX unsigned
#define TPP_MAX_FREQ_NUM 600
#define TPP_MAX_BONDS 4

#ifdef TPP_UNIT_TEST

 #define BOOST_TEST_BUILD_INFO yes
 #define BOOST_TEST_REPORT_LEVEL detailed
 #define BOOST_TEST_REPORT_FORMAT HRF
 #define BOOST_TEST_LOG_LEVEL all

    #include <boost/test/unit_test.hpp>
    #include <boost/test/unit_test_monitor.hpp>
    #include <boost/bind.hpp>
    #include <boost/test/framework.hpp>
    #include <boost/test/unit_test_log.hpp>
    #include <boost/test/results_collector.hpp>
    #include <boost/test/results_reporter.hpp>
    #include <boost/test/test_tools.hpp>
       using boost::unit_test::test_suite;
       using boost::bind;
       using boost::ref;
       namespace ut = boost::unit_test; 
 #define _BLOB_ \
                perror("UNIT_TEST is turned on so you should make all the BLOBS..\n"); \
                assert(0);
       
#else
 #include <cassert>
 #define _BLOB_(x) cout << "[ UNDER CONSTRUCTION - " << #x << " ]" << endl;
 // overloading of BOOST_TEST operators
 #define BOOST_CHECK(x) if (! ( x ) ) \
      printf("CHECK failed at procedure in %s at line %d\n --> %s", \
              __FILE__, __LINE__, #x);
 #define BOOST_REQUIRE(x) assert ( x );
 #define BOOST_FAIL(x) perror(#x); abort();
 #define BOOST_ERROR(x) perror(#x); abort();
 
 #ifdef ALLOW_WARNINGS
  #define BOOST_WARN(x) perror(#x);
 #else
  #define BOOST_WARN(x) ((void)0)
 #endif

#endif

namespace tpp {

using namespace boost::multi_index;
using namespace boost::posix_time;
using namespace boost::numeric;
using boost::array;
using boost::format;
using boost::numeric_cast;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::multi_index_container;
using OpenBabel::OBMol;

// atom and atom_array definition
// typedef struct { double x,y,z; } tpp_point;
typedef ublas::bounded_vector<double,3> t_point;

struct t_atom {
  TPP_INDEX                index; // place in vector array
  TPP_INDEX             oldindex;
  string               atom_type; 
  string              atom_type2;
  string               atom_name;
  string               old_aname;
  unsigned char          ncharge; // nuclear charge
  string                res_name;
  string                 comment;
  unsigned char           mol_id;
  t_point                  coord; 
  array<TPP_INDEX,TPP_MAX_BONDS>  connects;
  unsigned char     num_connects;
  double                  charge;
  double                    mass;
  TPP_INDEX                c_gnr;
  string                  qmname;
};

typedef multi_index_container<
         t_atom,
         indexed_by<
          ordered_unique<member<t_atom, TPP_INDEX, &t_atom::index> >, // key index (like array)
          ordered_non_unique<member<t_atom, string, &t_atom::atom_name> >, // key by name
          ordered_non_unique<member<t_atom, TPP_INDEX, &t_atom::c_gnr> > // key by charge group
         >
        > t_atom_array;

// I-O format
typedef enum { TPP_IF_PDB=0, TPP_IF_GRO=1, TPP_IF_G96=2, 
               TPP_IF_GAMOPT=3, TPP_IF_GAMHESS=4, TPP_IF_GAMSP=5 } 
        t_iformat;
typedef enum { TPP_OF_PDB=0, TPP_OF_GRO=2, TPP_OF_G96=3,
              TPP_OF_GAMIN=4 } t_oformat;

// internal coordinates
typedef enum {TPP_IC_BON=1, TPP_IC_ANG=2, TPP_IC_DIH=3, TPP_IC_CHIR=4} t_int_coord_type;

typedef struct { 
  t_int_coord_type type;
  TPP_INDEX i,j,k,l;
  string defname;
} t_int_coord;

typedef vector<t_int_coord> t_internals_array;

// topology parameters definition
typedef enum { TPP_TTYPE_BON=1, TPP_TTYPE_ANG=2, TPP_TTYPE_RBDIH=3, 
               TPP_TTYPE_IMPDIH=4,TPP_TTYPE_SYMDIH=5, TPP_TTYPE_PAIR=6,
               TPP_TTYPE_EXCL=7, TPP_TTYPE_SPECIMP=8 } t_ttype_type;

struct t_top_coord {
  string defname;
  t_ttype_type  type;
  short f; // f = -1 (means parameter lack)
  double c0,c1,c2,c3,c4,c5;/*
  t_top_coord(): defname(""), type(1),f(-1),c0(0),c1(0),c2(0),c3(0),c4(0),c5(0) {;}
  t_top_coord(const t_top_coord &_t): defname(_t.defname), type(_t.type),f(_t.f),c0(_t.c0),c1(_t.c1),c2(_t.c2),c3(_t.c3),c4(_t.c4),c5(_t.c5) {;}  */
};

struct t_top_element {
  string defname;
  TPP_INDEX i,j,k,l;
};

typedef multi_index_container< /* DEFINITIONS CONTAINER */
           t_top_coord, 
           indexed_by<
               ordered_unique<member<t_top_coord, string, &t_top_coord::defname> >,       // key by defname
               ordered_non_unique<member<t_top_coord, t_ttype_type, &t_top_coord::type> >,// key by directive
               ordered_non_unique<member<t_top_coord, short, &t_top_coord::f> >           // key by function
               > 
        > t_top_map;

typedef multi_index_container< /* TOPOLOGY CONTAINER */
           t_top_element,
           indexed_by<
            sequenced<>, // array key
            ordered_non_unique<member<t_top_element, string, &t_top_element::defname> >
           >             // key by defname
        > t_top_array;

// internal topology definition
class t_topology {
  public:
  t_top_map parameters;
  t_top_array elements;
  t_atom_array atoms;
  OBMol          mol;

  unsigned short nrexcl;
  string     ffinclude;
  string      res_name;
  string          name;

  /* serialization support */

  friend class boost::serialization::access;
    
  template<class Archive>
  void serialize(Archive& ar,const unsigned int)
  {
    ar & BOOST_SERIALIZATION_NVP(parameters);
    ar & BOOST_SERIALIZATION_NVP(elements);
    ar & BOOST_SERIALIZATION_NVP(atoms);
    ar & BOOST_SERIALIZATION_NVP(nrexcl);
    ar & BOOST_SERIALIZATION_NVP(ffinclude);
    ar & BOOST_SERIALIZATION_NVP(res_name);
  }
};

// optimize constant list
typedef vector<double*>  t_optimize_list;

// parameters for everything
typedef map<string, string> t_input_params;
typedef pair<string,string> t_input_param;
// fast work with params
template<typename T>
bool PARAM_EXISTS(const t_input_params &pars, const T x) {
  return (pars.find(x) != pars.end());
}
template<typename T>
const string PARAM_READ(const t_input_params &pars, const T x) {
  return (pars.find(x) != pars.end()) ? ((pars.find(x))->second) : string("");
}
template<typename T1, typename T2>
void PARAM_ADD(t_input_params &pars, const T1 par, const T2 val) {
   BOOST_CHECK( (pars.insert(t_input_param(par,val))).second );
}
template<typename T>
void PARAM_DEL(t_input_params &pars, const T par) {
   BOOST_CHECK(PARAM_EXISTS(pars,par)); 
   pars.erase( pars.find(par) );
}

// -- for example for command line :)
extern t_input_params cmdline;

class t_runtime {
  private:
    FILE *cash;
    FILE *log;
  public:
    t_runtime(const char *, const char *); // opening files log & cash
    void log_write(const char *); // string writing to log
    void log_write(string s) {log_write(s.c_str());}
    void cash_write(const char *, unsigned); // binary writing to cash
    ~t_runtime();
};
extern t_runtime runtime;

class t_exception {
  protected:
   string mesg;
   t_input_params pars;
  public:
   t_exception(): mesg("Undefined exception.") {;}
   t_exception(const char *, t_input_params &);
   t_exception(const char *s): mesg(s) {;}
   virtual string operator [] (const char *s) const { return string(PARAM_READ(pars,s)); }
   virtual string operator [] (const string &s) const { return string(PARAM_READ(pars,s)); }

   // rebuild all the files if you change pure-virtual function
   virtual void fix_log() const {
     std::ostringstream os;
     os << "TPP catched exception!\n";
     os << format("***** from %1% -> %2%\n") % PARAM_READ(pars, "classname") % PARAM_READ(pars, "procname");
     os << format("***** ===[ %1% ]===\n") % PARAM_READ(pars, "error");
     if (PARAM_EXISTS(pars, "line"))
         os << "***** parsing line: #" << PARAM_READ(pars, "line") << endl;
     os << "***** " << mesg << endl;
     runtime.log_write(os.str()); 
   }
};

// TODO: rewrite via MYSQLPP::QUERY interface
class t_db_exception: public t_exception {
  public:
   t_db_exception(const char *a, t_input_params &b): t_exception(a,b) {}
   virtual void fix_log() const {
     std::ostringstream os;
     os << "TPP catched exception!\n";
     os << format("***** from %1% -> %2%\n") % PARAM_READ(pars, "classname") % PARAM_READ(pars, "procname");
     os << format("***** ===[ %1% ]===\n") % PARAM_READ(pars, "error");
     if (PARAM_EXISTS(pars, "line"))
         os << "***** parsing line: #" << PARAM_READ(pars, "line") << endl;
     os << "***** SQL error: #" << PARAM_READ(pars, "sql_error") << endl;
     os << "***** SQL query: #" << PARAM_READ(pars, "sql_query") << endl;
     os << "***** " << mesg << endl;
     runtime.log_write(os.str()); 

   }
};

} // end namespace tpp

template<typename T>
ostream& operator<< (ostream& out, const vector<T>& v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last) 
            out << ", ";
    }
    out << "]";
    return out;
}


#define MYSQLPP_RESULT mysqlpp::StoreQueryResult
//#define MYSQLPP_RESULT mysqlpp::Result

#endif
