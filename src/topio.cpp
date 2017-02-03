#include "topio.hpp"
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

// anonymous for internal constants
namespace {
// comment for writing to the .itp/.rtp head
const char * top_comment =
  "; ----------------------------------------------\n"
  "; TPP - topology generator version " PACKAGE_VERSION " \n"
  "; created by Erg Research Group\n"
  "; MSU, Biology Faculty, Department of Biophysics\n"
  "; ----------------------------------------------\n"
  "; ATTENTION! Do not forget to use the proper version\n"
  "; of the force field fork (not less than revision). \n"
  "; Watch for corresponding force field at: \n"
  ";            bitbucket.com/comcon1\n"
  "; ----------------------------------------------\n"
  "; Please ascertain that the topology is valid. We \n"
  "; do not guarantee that. If you find that something\n"
  "; is wrong, please report us to " PACKAGE_BUGREPORT "\n";
}

namespace tpp {

  using OpenBabel::OBConversion;
  using OpenBabel::OBMol;
  using OpenBabel::OBMolAtomIter;
  using OpenBabel::OBAtomAtomIter;
  using namespace phoenix;

  void save_topology_rtp(const Topology &tp, const char *fname) {
    // test if file exists
    TPPD << format("Trying to write RTP-topology into '%s'.")  % fname;
    fstream out(fname, ios::out);
    if (!out.is_open()) {
      TPPE << "Fail to open file for write.";
      Parameters params;
      params.add("procname", "tpp::save_topology_rtp");
      params.add("error", "invalid filename");
      params.add("filename", fname);
      Exception e("Can't open specified file for write.", params);
      throw e;
    }
    // header
    out << format("\
  %1%\
  [ %2$4s ]\n\
  ") % top_comment % tp.res_name;

    // atoms specification
    out << "\n[ atoms ]\n";
    for (auto &it: tp.atoms) {
      out << format("  %1$-4s   %2$-4s   %3$6.3f   %4$3d \n") %
             it.atom_name.c_str() % it.atom_type.c_str() % it.charge %
             (int)(it.c_gnr);
    }
    out << flush;

    // bonds
    if (tp.parameters.get<1>().find(TPP_TTYPE_BON) != tp.parameters.get<1>().end()) {
      out << "\n[ bonds ]\n";
      for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
           it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0)
          if (it->f == -1) {
            out << format(" %1$-4s   %2$-4s   bndTPP_%3$s_%4$s \n")
                    % tp.atoms.find((int)it0->i)->atom_name.c_str()
                    % tp.atoms.find((int)it0->j)->atom_name.c_str()
                    % tp.atoms.find((int)it0->i)->atom_type2.c_str()
                    % tp.atoms.find((int)it0->j)->atom_type2.c_str();
          } else {
            out << format(" %1$-4s   %2$-4s \n")
                    % tp.atoms.find((int)it0->i)->atom_name.c_str()
                    % tp.atoms.find((int)it0->j)->atom_name.c_str();
          }
    }

    // angles
    bool flag = true;
    for (auto it = tp.parameters.get<2>().lower_bound(-1);
        it != tp.parameters.get<2>().upper_bound(-1); ++it)
      if (it->type == TPP_TTYPE_ANG)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0) {
          if (flag) {
            out << "\n[ angles ]\n";
            flag = false;
          }
          out << format("%1$4s %2$4s %3$4s dhTPP_%4$s_%5$s_%6$s \n")
                  % tp.atoms.find((int)it0->i)->atom_name.c_str()
                  % tp.atoms.find((int)it0->j)->atom_name.c_str()
                  % tp.atoms.find((int)it0->k)->atom_name.c_str()
                  % tp.atoms.find((int)it0->i)->atom_type2.c_str()
                  % tp.atoms.find((int)it0->j)->atom_type2.c_str()
                  % tp.atoms.find((int)it0->k)->atom_type2.c_str();
        }

    // dihedrals
    flag = true;
    for (auto it = tp.parameters.get<2>().lower_bound(-1);
        it != tp.parameters.get<2>().upper_bound(-1); ++it)
      if (it->type == TPP_TTYPE_RBDIH)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0) {
          if (flag) {
            out << "\n[ dihedrals ]\n";
            flag = false;
          }
          out << format("%1$4s %2$4s %3$4s %4$4s dhTPP_%5$s_%6$s_%7$s_%8$s \n")
                  % tp.atoms.find((int)it0->i)->atom_name.c_str()
                  % tp.atoms.find((int)it0->j)->atom_name.c_str()
                  % tp.atoms.find((int)it0->k)->atom_name.c_str()
                  % tp.atoms.find((int)it0->l)->atom_name.c_str()
                  % tp.atoms.find((int)it0->i)->atom_type2.c_str()
                  % tp.atoms.find((int)it0->j)->atom_type2.c_str()
                  % tp.atoms.find((int)it0->k)->atom_type2.c_str()
                  % tp.atoms.find((int)it0->l)->atom_type2.c_str();
        }

    // impropers
    flag = true;
    for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_SPECIMP);
        it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SPECIMP); ++it)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
          it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0) {
          if (flag) {
            out << "\n[ impropers ]\n";
            flag = false;
          }
           out << format("%1$4s %2$4s %3$4s %4$4s %5$s \n")
                  % tp.atoms.find((int)it0->i)->atom_name.c_str()
                  % tp.atoms.find((int)it0->j)->atom_name.c_str()
                  % tp.atoms.find((int)it0->k)->atom_name.c_str()
                  % tp.atoms.find((int)it0->l)->atom_name.c_str()
                  % it->defname.c_str();
        }
    out << "; end of TPPMKTOP topology" << endl;
    out.close();
    TPPD << "Topology has been succesfully written.";
  } // save_topology_rtp


  void save_topology(const Topology &tp, const char *fname, bool ncf) {
    TPPD << format("Trying to write ITP-topology into '%s'.")  % fname;
    fstream out(fname, ios::out);
    if (!out.is_open()) {
      TPPE << "Fail to open file for write.";
      Parameters params;
      params.add("procname", "tpp::save_topology");
      params.add("error", "invalid filename");
      params.add("filename", fname);
      Exception e("Can't open specified file for write.", params);
      throw e;
    }
    // header
    out << format("\
  %1%\
  ; ----------------------------------------------\n\
  ; Topology was prepared for use with the force field:\n\
  ;         %2%\n\
  \n\
  [ moleculetype ]\n\
   %3$4s %4$1d\n\
  \n\
  ; Force constant parameters\n") % top_comment % tp.ffinfo % tp.res_name % (int)(tp.nrexcl);
    // force constants parameters '#define's
    for (auto it = tp.parameters.get<1>().begin();
          it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it) {
      if ( (it->f != -1) ^  ncf ) {
        switch (it->type) {
         case TPP_TTYPE_BON:   out << format("#define %1$-25s %2$1d %3$8.3f %4$9.2e ;1\n")
                                % it->defname % it->f % it->c0 % it->c1; break;
         case TPP_TTYPE_ANG:   out << format("#define %1$-25s %2$1d %3$6.1f %4$8.1f    ;2\n")
                               % it->defname % it->f % it->c0 % it->c1; break;
         case TPP_TTYPE_RBDIH: out << format("#define %1$-25s 3 %2$+5.1f %3$+5.1f %4$+5.1f %5$+5.1f %6$+5.1f %7$+5.1f ;3\n")
                               % it->defname % it->c0 % it->c1 % it->c2 % it->c3 % it->c4 % it->c5; break;
         case TPP_TTYPE_IMPDIH:out << format("#define %1$-25s 2 %2$5.1f %3$6.1f ;4\n")
                               % it->defname % it->c0 % it->c1; break;
         case TPP_TTYPE_SYMDIH:out << format("#define %1$-25s 1 %2$5.1f %3$6.1f %4$1d ;5\n")
                               % it->defname % it->c0 % it->c1 % numeric_cast<unsigned>(it->c2);
        };
      }

    }
    // atoms specification
    out << "\n[ atoms ]\n";
    for (auto &it: tp.atoms) {
      out << format("%1$3d  %2$-10s 1  %3$-4s  %4$-4s  %5$2d  %6$+6.3f  %7$10.6f ; [%8$3s] %9$s\n") %
             (int)(it.index) % it.atom_type.c_str() % tp.res_name.c_str() % it.atom_name.c_str() %
             (int)(it.c_gnr) % it.charge % it.mass % it.atom_type2.c_str() % it.comment;
    }
    out << flush;
    // bonds
    if (tp.parameters.get<1>().find(TPP_TTYPE_BON) != tp.parameters.get<1>().end()) {
      out << "\n[ bonds ]\n";
      for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
           it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0)
          if (ncf && (it->f != -1) )
            out << format("%1$3d %2$3d  %3$2d\n") % (int)it0->i % (int)it0->j % it->f;
          else
            out << format("%1$3d %2$3d  %3$-15s\n") % (int)it0->i % (int)it0->j % it0->defname;
    }
    // angles
    if (tp.parameters.get<1>().find(TPP_TTYPE_ANG) != tp.parameters.get<1>().end()) {
      out << "\n[ angles ]\n";
      for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
           it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0)
          if (ncf && (it->f != -1) )
            out << format("%1$3d %2$3d %3$3d  %4$2d\n") % (int)it0->i % (int)it0->j % (int)it0->k % it->f;
          else
            out << format("%1$3d %2$3d %3$3d  %4$-15s\n") % (int)it0->i % (int)it0->j % (int)it0->k % it0->defname;
    }
    // dihedrals
    if ( ( tp.parameters.get<1>().find(TPP_TTYPE_RBDIH) != tp.parameters.get<1>().end() ) ||
        //TODO: incorporate IMPDIH here?
         ( tp.parameters.get<1>().find(TPP_TTYPE_SYMDIH) != tp.parameters.get<1>().end() )
        ) {
      out << "\n[ dihedrals ]\n";
      for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
           it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0)
          if (ncf && (it->f != -1) )
            out << format("%1$3d %2$3d %3$3d %4$3d  %5$2d\n") % (int)it0->i % (int)it0->j % (int)it0->k % (int)it0->l % it->f;
          else
            out << format("%1$3d %2$3d %3$3d %4$3d  %5$-15s\n") % (int)it0->i % (int)it0->j % (int)it0->k % (int)it0->l % it0->defname;
    }
    // impropers
    if ( tp.parameters.get<1>().find(TPP_TTYPE_SPECIMP) != tp.parameters.get<1>().end() ) {
      out << "\n[ dihedrals ]\n";
      for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_SPECIMP);
           it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SPECIMP); ++it)
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0)
          // impropers are always with defines (!)
            out << format("%1$3d %2$3d %3$3d %4$3d  %5$-15s\n") % (int)it0->i % (int)it0->j % (int)it0->k % (int)it0->l % it0->defname;
          //TODO: default improper type should be incorporated
    }
    // pairs
    if ( tp.parameters.get<1>().count(TPP_TTYPE_PAIR) > 0 ) {
      out << "\n[ pairs ]\n";
      for (auto it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_PAIR);
           it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_PAIR); ++it) {
        out << ";\n";
        for (auto it0 = tp.elements.get<1>().lower_bound(it->defname);
             it0 != tp.elements.get<1>().upper_bound(it->defname); ++it0) {
          if (it0->defname == "ONE_PAIR")
            out << format("%1$3d %2$3d  %3$1d\n") % (int)it0->i % (int)it0->j % it->f;
          else
            out << format("%1$3d %2$3d  %3$1d %4$10.5f %5$10.5f\n") % (int)it0->i % (int)it0->j % it->f % it->c0 % it->c1;
        }
      }
    }
    // exclusion
    if ( tp.parameters.get<1>().find(TPP_TTYPE_EXCL) != tp.parameters.get<1>().end() ) {
      out << "\n[ exclusions ]\n";
      string defname = (tp.parameters.get<1>().find(TPP_TTYPE_EXCL))->defname; // in map there is only one exclusion
        for (auto it0 = tp.elements.get<1>().lower_bound(defname);
             it0 != tp.elements.get<1>().upper_bound(defname); ++it0)
          out << format("%1$3d %2$3d\n") % it0->i % it0->j;
    }
    out << "; end of TPPMKTOP topology" << endl;
    out.close();
    TPPD << "Topology has been sucessfully written.";
  } // save_topology

  // save information about lacking topology
  void save_lack(const Topology &tp, const char *fname) {
        TPPI << format("TPP will write %1$d lack parameters to %2$s.\n")
          % tp.parameters.get<2>().count(-1) % fname;
        std::ofstream qalcfile(fname, ios::out);
        if (!qalcfile.is_open()) {
          TPPE << "Fail to open file for write.";
          Parameters params;
          params.add("procname", "tpp::save_lack");
          params.add("error", "invalid filename");
          params.add("filename", fname);
          Exception e("Can't open specified file for write.", params);
          throw e;
        }
        qalcfile << "; TPP topology lack\n";
        for (auto typit = tp.parameters.get<2>().lower_bound(-1);
           typit != tp.parameters.get<2>().upper_bound(-1); ++typit) {
          switch (typit->type) {
            case TPP_TTYPE_BON:
             qalcfile << format("#define %1$s 1 0.00 0.00 ;1\n") % typit->defname;
             break;
            case TPP_TTYPE_ANG:
             qalcfile << format("#define %1$s 1 0.00 0.00 0.00 ;2\n") % typit->defname;
             break;
            case TPP_TTYPE_IMPDIH:
             qalcfile << format("#define %1$s 2 0.00 0.00 0.00 ;3\n") % typit->defname;
             break;
            case TPP_TTYPE_RBDIH:
             qalcfile << "; May be overwritten with pairs \n";
             qalcfile << format("#define %1$s 3 0.00 0.00  0.00 0.00 0.00 ;4\n") % typit->defname;
             break;
          };
        }
        qalcfile << flush;
        qalcfile.close();
        TPPD << "Lack file has been written.";
  } // end save_lack

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
          Parameters params;
          params.add("procname", "tpp::load_lack");
          params.add("error", "invalid filename");
          params.add("filename", fname);
          Exception e("Can't open specified file for read.", params);
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
         Parameters pars0;
         pars0.add("error", "Parsing error");
         throw Exception("Topology lack parsing error!", pars0);
   }
   return;

  } // end load_lack



  void load_topology(Topology &tp, const char *fname) {
    TPPD << format("Trying to read ITP file into topology from '%s'.") % fname;
    { // checking file readability
      fstream inf(fname, ios::in);
      if (!inf.is_open()) {
          TPPE << "Fail to open file for read.";
          Parameters params;
          params.add("procname", "tpp::load_topology");
          params.add("error", "invalid filename");
          params.add("filename", fname);
          Exception e("Can't open specified file for read.", params);
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
         Parameters pars0;
         pars0.add("error", "Parsing error");
         throw Exception("Topology parsing error!", pars0);
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



// also serialization should be included

} // end namespace

