#include "topreader.hpp"

#include "exceptions.hpp"
#include <cassert>
#include "logger.hpp"

#include "strutil.hpp"

#include <boost/format.hpp>

#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <cctype>
#include <iterator>

#define PHOENIX_LIMIT 5

#if HAVE_BOOST_SPIRIT_CORE_HPP
	#include <boost/spirit/attribute.hpp>
	#include <boost/spirit/phoenix/primitives.hpp>
	#include <boost/spirit/phoenix/operators.hpp>
	#include <boost/spirit/phoenix/functions.hpp>
	#include <boost/spirit/phoenix/binders.hpp>
	#include <boost/spirit/phoenix/casts.hpp>
	#include <boost/spirit/utility/regex.hpp>

#elif HAVE_BOOST_SPIRIT_HOME_CLASSIC_CORE_HPP
	#include <boost/spirit/home/classic/attribute.hpp>
	#include <boost/spirit/home/classic/phoenix/primitives.hpp>
	#include <boost/spirit/home/classic/phoenix/operators.hpp>
	#include <boost/spirit/home/classic/phoenix/functions.hpp>
	#include <boost/spirit/home/classic/phoenix/binders.hpp>
	#include <boost/spirit/home/classic/phoenix/casts.hpp>
	#include <boost/spirit/home/classic/utility/regex.hpp>
#endif

#include "lexical.hpp"

#if ENABLE_GAMESS_FEATURES
	#include <boost/numeric/ublas/matrix_proxy.hpp>
#endif

#define OPENBABEL
#ifdef DEBUG
    #define BOOST_SPIRIT_DEBUG
#endif // DEBUG

namespace tpp {


  using std::string;
  using std::pair;
  using std::cout;
  using std::endl;
  using std::ios;
  using std::flush;
  using std::ostream;
  using std::fstream;
  using std::ostringstream;

  using boost::numeric_cast;
  using boost::format;

  #if ENABLE_GAMESS_FEATURES
  typedef ublas::matrix<double> umd;
  static member_function_ptr<umd::reference, umd,  umd::size_type, umd::size_type, umd::const_reference>
   ins_elem = &umd::insert_element;

  void load_hessian(ublas::matrix<double>& mtx, const char *fname) {
     {
    TPPD << format("Trying to read %s into hessian matrix.") % fname;
    fstream inf(fname, ios::in);
    if (!inf.is_open()) {
      TPPD << "Fail to open file for read.";
      BOOST_CHECK(0);
    }
    inf.close();
     }

    lex::iterator first(fname);
    lex::iterator last = first.make_end();
    lex::rule_r soms, define, common;
    int num_x = 0, num_y = 0;
    int x = 1, y = 1;
    soms = ( (*(anychar_p - eol_p)) - ( +blank_p >> str_p("CARTESIAN FORCE CONSTANT MATRIX") >> *(anychar_p - eol_p) ) )  >> eol_p;
    common = +soms >>
      +blank_p >> str_p("CARTESIAN FORCE CONSTANT MATRIX")[var(x) = 0][var(y) = 0]>> *(anychar_p - eol_p) >> eol_p
      >> *(anychar_p - eol_p) >> eol_p >>
      // starting one block
      +( *blank_p >> eol_p >>
        (+blank_p >> uint_p[assign_a(num_x)] >> !(+blank_p >> uint_p) >> *blank_p >> eol_p)[&lex::echo::val_print] >>
        (+blank_p >> +alnum_p >> !(+blank_p >> +alnum_p) >> *blank_p >> eol_p)[&lex::echo::val_print] >>
        *(anychar_p - eol_p) >> eol_p >> // " x y z x y z" string
        repeat_p(3,boost::spirit::classic::more)
        [
        +blank_p >> !( uint_p[assign_a(num_y)][var(y) = var(num_y)*3-3] >> +blank_p >> +alnum_p >> +blank_p )[&lex::echo::val_print]
                 >> (ch_p('X') ^ 'Y' ^ 'Z')[increment_a(y)][var(x) = var(num_x)*3-3]
                 >> (repeat_p(3,6)[ !blank_p >> strict_real_p[increment_a(x)][ins_elem(var(mtx),var(x)-1,var(y)-1, arg1)] ])
                 >> eol_p
        ]
      )
      // ending block
      >> +soms;
   parse_info<lex::iterator> info;
   //BOOST_SPIRIT_DEBUG_NODE( common );
   try {
     info = parse(first,last, common);
   } catch (parser_error_base &e) {
        cout << "Parse error. Check file format (GAMESS OUTPUT, RUNTYP=HESSIAN).\n";
   }
   cout << format("%1$4d %2$4d %3$4d %4$4d\n") % x % y % num_x % num_y;
   cout << "Current hesse matrix is: \n";

   if (!info.full) {
        cout << format("\
  Error in parsing file '%1%' catched:  \n\
  ---------------------------------\n") % fname;
          lex::iterator err_it(first);
          std::advance(err_it, info.length);
    std::copy(first, err_it,
                   std::ostream_iterator<lex::iterator::value_type>(cout)
                  );
          cout << "*\n\
  ---------------------------------\n\
  * - the point of error\n";
         Parameters pars0;
         PARAM_ADD(pars0, "error", "Parsing error");
         throw Exception("GAMESS output parsing error!", pars0);
   }

   // symmetry matrix ;))
   for (int i=0; i<mtx.size1(); ++i)
    for (int j=0; j<i; ++j)
     mtx(i,j) = mtx(j,i);
   //
   ostringstream os;
   os << "Matrix hessian from GAMESS-file:\n";
   for (int i=0; i<mtx.size1(); ++i)
     os << ublas::row(mtx,i) << endl;
   TPPD << os.str();
   return;
  }
  #endif // CONDITIONAL COMPILE: ENABLE_GAMESS_FEATURES

  // loading topology parameters from lack-file
  void load_lack(Topology &tp, const char *fname) {
    TPPD << format("Trying to read lack-file into topology from '%s'.") % fname;
    { // checking file readability
      fstream inf(fname, ios::in);
      if (!inf.is_open()) {
          TPPE << "Fail to open file for read.";
          Exception e("Can't open specified file for read.");
          e.add("procname", "tpp::load_lack");
          e.add("error", "invalid filename");
          e.add("filename", fname);
          throw e;
      }
      inf.close();
     }
    lex::iterator first(fname);
    lex::iterator last = first.make_end();
    lex::rule_r clrs, coms, define, common;
    TopCoord _tc;

    coms = *blank_p >> ch_p(';') >> (*(anychar_p - eol_p))[&lex::echo::val_print] >> eol_p;
    clrs = *blank_p >> eol_p;
    //define string
    define  = *blank_p >> str_p("#define") >> +blank_p
                       >> (+(alnum_p | '_'))[assign_a(_tc.defname)][&lex::echo::val_print] >> +blank_p
                       >> uint_p[assign_a(_tc.f,-1)] >> +blank_p
                       >> real_p[assign_a(_tc.c0)] >> +blank_p
                       >> real_p[assign_a(_tc.c1)] >> +blank_p
                       >> !(real_p[assign_a(_tc.c2)] >> +blank_p
                         >>   !( real_p[assign_a(_tc.c3)] >> +blank_p
                              >> real_p[assign_a(_tc.c4)] >> +blank_p
                              >> real_p[assign_a(_tc.c5)] >> +blank_p
                             )
                           ) >> ch_p(';') >> uint_p[lex::assign_cast_a(_tc.type)] >> *blank_p
                    >> eol_p[lex::insert_mi_a(tp.parameters,_tc)];
   common = +(clrs ^ coms ^ define);

   parse_info<lex::iterator> info;

   BOOST_SPIRIT_DEBUG_NODE( common );
   try {
     info = parse(first,last, common);
   } catch (const parser_error_base &e) {
        cout << "Parse error. Check file format.\n";
   }
   cout << format("\
  %1$4d define-parameters found\n") % tp.parameters.size() << endl;
   if (!info.full) {
        cout << format("\
  Error in parsing file '%1%' catched:  \n\
  ---------------------------------\n") % fname;
          lex::iterator err_it(first);
          std::advance(err_it, info.length);
    std::copy(first, err_it,
                   std::ostream_iterator<lex::iterator::value_type>(cout)
                  );
          cout << "*\n\
  ---------------------------------\n\
  * - the point of error\n";
          cout << "Parsed parameters: " << endl;
          cout << lex::outs.str() << endl;
         Exception e("Topology lack parsing error!");
         e.add("error", "Parsing error");
         throw e;
   }
   return;

  } // end load_lack



  void load_topology(Topology &tp, const char *fname) {
    TPPD << format("Trying to read ITP file into topology from '%s'.") % fname;
    { // checking file readability
      fstream inf(fname, ios::in);
      if (!inf.is_open()) {
          TPPE << "Fail to open file for read.";
          Exception e("Can't open specified file for read.");
          e.add("procname", "tpp::load_topology");
          e.add("error", "invalid filename");
          e.add("filename", fname);
          throw e;
      }
      inf.close();
     }
    // header
    lex::iterator first(fname);
    lex::iterator last = first.make_end();
    lex::rule_r clrs, coms, define, atom, bond, angl, dihe, header,
      pair, atoms, moltype, bonds, angles, dihedrals, pairs, common;
    TopCoord _tc;
    Atom _at;
    TopElement _te;

    coms = *blank_p >> ch_p(';') >> (*(anychar_p - eol_p))[&lex::echo::val_print] >> eol_p;
    clrs = *blank_p >> eol_p;
    //define string
    define  = *blank_p >> str_p("#define") >> +blank_p
                       >> (+(alnum_p | '_'))[assign_a(_tc.defname)][&lex::echo::val_print] >> +blank_p
                       >> uint_p[assign_a(_tc.f)] >> +blank_p
                       >> real_p[assign_a(_tc.c0)] >> +blank_p
                       >> real_p[assign_a(_tc.c1)] >> +blank_p
                       >> !(real_p[assign_a(_tc.c2)] >> +blank_p
                         >>   !( real_p[assign_a(_tc.c3)] >> +blank_p
                              >> real_p[assign_a(_tc.c4)] >> +blank_p
                              >> real_p[assign_a(_tc.c5)] >> +blank_p
                             )
                           ) >> ch_p(';') >> uint_p[lex::assign_cast_a(_tc.type)] >> *blank_p
                    >> eol_p[lex::insert_mi_a(tp.parameters,_tc)];
    // atom string
    // if atoms already loaded from structure file
    if (!tp.atoms.size())
    atom = *blank_p >> uint_p[assign_a(_at.index)] >> +blank_p // atom index
                    >> (+(alnum_p | '_' | '='))[assign_a(_at.atom_type)] >> +blank_p // atom type OPLS
                    >> uint_p >> +blank_p // chain
                    >> (+alnum_p) >> +blank_p // residue
                    >> (+alnum_p)[assign_a(_at.atom_name)] >> +blank_p
                    >> uint_p[assign_a(_at.c_gnr)] >> +blank_p
                    >> real_p[assign_a(_at.charge)] >> +blank_p
                    >> real_p[assign_a(_at.mass)] >> +blank_p
                    >> ch_p(';') >> (+(anychar_p - eol_p))[assign_a(_at.comment)]
                  >> eol_p[lex::insert_mi_a(tp.atoms,_at)];
    else // if no atoms present in topology now
    atom = *blank_p >> uint_p[assign_a(_at.index)] >> +blank_p // atom index
                    >> (+(alnum_p | '_' | '='))[assign_a(_at.atom_type)] >> +blank_p // atom type OPLS
                    >> uint_p >> +blank_p // chain
                    >> (+alnum_p) >> +blank_p // residue
                    >> (+alnum_p)[assign_a(_at.atom_name)] >> +blank_p
                    >> uint_p[assign_a(_at.c_gnr)] >> +blank_p
                    >> real_p[assign_a(_at.charge)] >> +blank_p
                    >> real_p[assign_a(_at.mass)] >> +blank_p
                    >> ch_p(';') >> (+(anychar_p - eol_p))[assign_a(_at.comment)]
                  >> eol_p[lex::update_mi_a(tp.atoms,_at)];
    // bond string
    bond = *blank_p >> uint_p[assign_a(_te.i)] >> +blank_p
                    >> uint_p[assign_a(_te.j)] >> +blank_p
                    >> (+(alnum_p | '_'))[assign_a(_te.defname)] >> *blank_p
                >> eol_p[lex::insert_mi_a(tp.elements.get<1>(), _te)];
    // angle string
    angl = *blank_p >> uint_p[assign_a(_te.i)] >> +blank_p
                    >> uint_p[assign_a(_te.j)] >> +blank_p
                    >> uint_p[assign_a(_te.k)] >> +blank_p
                    >> (+(alnum_p | '_'))[assign_a(_te.defname)] >> *blank_p
                >> eol_p[lex::insert_mi_a(tp.elements.get<1>(), _te)];
    // dihedral string
    dihe = *blank_p >> uint_p[assign_a(_te.i)] >> +blank_p
                    >> uint_p[assign_a(_te.j)] >> +blank_p
                    >> uint_p[assign_a(_te.k)] >> +blank_p
                    >> uint_p[assign_a(_te.l)] >> +blank_p
                    >> (+(alnum_p | '_'))[assign_a(_te.defname)] >> *blank_p
                >> eol_p[lex::insert_mi_a(tp.elements.get<1>(), _te)];
    pair = *blank_p >> uint_p[assign_a(_te.i)] >> +blank_p
                    >> uint_p[assign_a(_te.j)] >> +blank_p
                    >> ch_p('1')[assign_a(_te.defname,"ONE_PAIR")] >> *blank_p
                    >> eol_p[lex::insert_mi_a(tp.elements.get<1>(), _te)];
    // moleculetype block
    moltype = *blank_p >> ch_p('[') >> *blank_p >> str_p("moleculetype") >> *blank_p >> ch_p(']') >> *blank_p >> eol_p
      >> *blank_p >> (+alnum_p)[assign_a(tp.res_name)] >> +blank_p >> uint_p[assign_a(tp.nrexcl)] >> *blank_p >> eol_p;
    atoms =  *blank_p >> ch_p('[') >> *blank_p >> str_p("atoms") >> *blank_p >> ch_p(']') >> *blank_p >> eol_p
          >> +atom;
    bonds =  *blank_p >> ch_p('[') >> *blank_p >> str_p("bonds") >> *blank_p >> ch_p(']') >> *blank_p >> eol_p
      >> +bond;
    angles = *blank_p >> ch_p('[') >> *blank_p >> str_p("angles") >> *blank_p >> ch_p(']') >> *blank_p >> eol_p
      >> +angl;
    dihedrals =  *blank_p >> ch_p('[') >> *blank_p >> str_p("dihedrals") >> *blank_p >> ch_p(']') >> *blank_p >> eol_p
      >> +dihe;
    pairs =  *blank_p >> ch_p('[') >> *blank_p >> str_p("pairs") >> *blank_p >> ch_p(']') >> *blank_p >> eol_p
      >> +pair;
    header = *blank_p >> str_p("#include") >> +blank_p >> ch_p('<')
      >> (+(anychar_p - '>' - eol_p))[assign_a(tp.ffinclude)][&lex::echo::val_print]
      >> '>' >> *blank_p >> eol_p;
   common = *(clrs ^ coms)
          >> +header
          >> *(clrs ^ coms)
          >> moltype
          >> *(clrs ^ coms)
          >> *(define)
          >> *(clrs ^ coms)
          >> atoms
          >> *(clrs ^ coms)
          >> !bonds
          >> *(clrs ^ coms)
          >> !angles
          >> *(clrs ^ coms)
          >> !dihedrals
          >> *(clrs ^ coms)
          >> !(pairs
                [assign_a(_tc.type, TPP_TTYPE_PAIR)]
                [assign_a(_tc.defname, "ONE_PAIR")]
                [assign_a(_tc.f, 1)]
                [lex::insert_mi_a(tp.parameters, _tc)]
               )
          >> *(clrs ^ coms);

   parse_info<lex::iterator> info;

   BOOST_SPIRIT_DEBUG_NODE( common );
   try {
     info = parse(first,last, common);
   } catch (const parser_error_base &e) {
        TPPE << "Parse error. Check file format.\n";
   }
   cout << format("\
  %1$4d define-parameters found\n\
  %2$4d topology elements found\n\
  %3$4d atoms found\n") % tp.parameters.size() % tp.elements.size() % tp.atoms.size() << endl;
   if (!info.full) {
        cout << format("\
  Error in parsing file '%1%' catched:  \n\
  ---------------------------------\n") % fname;
          lex::iterator err_it(first);
          std::advance(err_it, info.length);
          std::copy(first, err_it,
                   std::ostream_iterator<lex::iterator::value_type>(cout)
                  );
          cout << "*\n\
  ---------------------------------\n\
  * - the point of error\n";
          cout << "Parsed parameters: " << endl;
          cout << lex::outs.str() << endl;
         Exception e("Topology parsing error!");
         e.add("error", "Parsing error");
         throw e;
   }
   {
     ostringstream os;
     os << "Bonds readed: " << tp.parameters.get<1>().count(TPP_TTYPE_BON) << endl;
     os << "Angles readed: " << tp.parameters.get<1>().count(TPP_TTYPE_ANG) << endl;
     os << "RB dihedrals readed: " << tp.parameters.get<1>().count(TPP_TTYPE_RBDIH) << endl;
     os << "Pairs readed: " << tp.parameters.get<1>().count(TPP_TTYPE_PAIR) << endl;
     TPPD << os.str();
   }
   return;
  }

  void check_topology(Topology &tp) {

    // header

    ;
  }

}
