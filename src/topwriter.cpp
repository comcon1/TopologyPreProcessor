#include "topwriter.hpp"
#include "exceptions.hpp"
#include <cassert>
#include "logger.hpp"

#include "strutil.hpp"

#include <boost/format.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <cctype>
#include <iterator>

#define OPENBABEL

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

  using OpenBabel::OBConversion;
  using OpenBabel::OBMol;
  using OpenBabel::OBMolAtomIter;
  using OpenBabel::OBAtomAtomIter;

  // save RTP format
  void TopologyWriter::saveRTP(const Topology &tp, const char *fname) {
    // test if file exists
    TPPD << format("Trying to write RTP-topology into '%s'.")  % fname;
    fstream out(fname, ios::out);
    if (!out.is_open()) {
      TPPE << "Fail to open file for write.";
      Exception e("Can't open specified file for write.");
      e.add("procname", "tpp::save_topology_rtp");
      e.add("error", "invalid filename");
      e.add("filename", fname);
      throw e;
    }
    // header
    out << format("\
  %1%\
; RTP works properly only with FF: %3%\n\
  \
[ %2$4s ]\n\
  ") % top_comment % tp.res_name % tp.ffinfo;

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
  } // saveRTP


  // save topology to file in GMX
  void TopologyWriter::saveITP(const Topology &tp, const char *fname, bool ncf) {
    TPPD << format("Trying to write ITP-topology into '%s'.")  % fname;
    fstream out(fname, ios::out);
    if (!out.is_open()) {
      TPPE << "Fail to open file for write.";
      Exception e("Can't open specified file for write.");
      e.add("procname", "tpp::save_topology");
      e.add("error", "invalid filename");
      e.add("filename", fname);
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
  void TopologyWriter::saveAbsentParametersITP(const Topology &tp, const char *fname) {
        TPPI << format("TPP will write %1$d lack parameters to %2$s.\n")
          % tp.parameters.get<2>().count(-1) % fname;
        std::ofstream qalcfile(fname, ios::out);
        if (!qalcfile.is_open()) {
          TPPE << "Fail to open file for write.";
          Exception e("Can't open specified file for write.");
          e.add("procname", "tpp::save_lack");
          e.add("error", "invalid filename");
          e.add("filename", fname);
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

// also serialization should be included

} // end namespace

