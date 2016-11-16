#ifndef TPP_GLOBAL_H
#define TPP_GLOBAL_H

#ifdef HAVE_CONFIG_H
#include "../config.h"
#else
#error !! YOU SHOULD RUN CONFIGURE SCRIPT !!
#endif


#include "paramset.hpp"

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
// spirit switch

#if HAVE_BOOST_SPIRIT_CORE_HPP
#define SPIRIT_HOME boost/spirit
#elif HAVE_BOOST_SPIRIT_HOME_CLASSIC_CORE_HPP
#define SPIRIT_HOME boost/spirit/home/classic
#endif

//
//	THIS SECTION MUST BE REMOVED
//
//	It pollutes global namespace!
//
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

//
//	END of removed section
//
//
#define TPP_GMX_EXTERNAL
#define TPP_MAX_ATOM_NUM 200
#define TPP_INDEX unsigned
#define TPP_MAX_FREQ_NUM 600
#define TPP_MAX_BONDS 4



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
  long int dbid;
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
  t_top_map  parameters;
  t_top_array  elements;
  t_atom_array    atoms;
  OBMol             mol;

  unsigned short nrexcl;
  string      ffinclude;
  string         ffinfo;
  string       res_name;
  string           name;

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
