/*! \file tppmktop.cpp
 *
 *	\brief This file provides executable for TPPMKTOP utility.
 *
 */

#include "global.hpp"
#include "core.hpp"
#include "exceptions.hpp"
#include "logger.hpp"
#include "topio.hpp"
#include "structio.hpp"
#include "db_base.hpp"
#include "tppnames.hpp"
#include "atom_definer.hpp"
#include "bond_definer.hpp"
#include "strutil.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/filesystem.hpp>

#include <boost/format.hpp>
#include <sstream>

namespace p_o = boost::program_options;
namespace bfs = boost::filesystem;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using boost::format;

void printHelp(p_o::options_description const&);
void printInfo();

double sumcharge(const tpp::Topology &);

int main(int argc, char * argv[]) {
  tpp::initiate_logging("tppmktop.log");
  string progname("Execution rules for TPPMKTOP ");
  progname = progname + VERSION;
  p_o::options_description desc(progname);
  p_o::variables_map vars;
  desc.add_options()
        ("input,i",
            p_o::value<std::string>()->required(),
            "Input filename (any format)")
        ("output,o",
            p_o::value<std::string>()->required(),
            "Output filename (itp format)")
        ("rtp-output,r",
            p_o::value<std::string>()->default_value(""),
            "Output filename (rtp format)")
        ("forcefield,f",
            p_o::value<std::string>()->required(),
            "Forcefield name")
        ("lack-file,l",
            p_o::value<std::string>()->default_value("lack.itp"),
            "Topology lack filename (default 'lack.itp')")
        ("sqlserver,s",
            p_o::value<std::string>()->default_value("localhost"),
            "Mysql-server adress (default 'localhost')")
        ("sqlport,t",
            p_o::value<unsigned>()->default_value(3306),
            "Server port (default '3306')")
        ("sqluser,u",
            p_o::value<std::string>()->default_value("tppuser"),
            "Database username (default 'tppuser')")
        ("sqlpassword,p",
            p_o::value<std::string>()->default_value("estatic"),
            "Mysql-password (default 'estatic')")
        ("nocalculate,n",
            p_o::value<bool>()->default_value(false)->implicit_value(false),
            "Create final topology (don't create lack-file)")
        ("max-bonds,m",
            p_o::value<bool>()->default_value(false)->implicit_value(false),
            "Maximize amount of bonds, angles and dihedrals by selecting other atom-types.")
        ("verbose,v",
            p_o::value<bool>()->default_value(false)->implicit_value(false),
            "Verbose mode")
      ("help,h", p_o::bool_switch()->default_value(false),
          "Print detailed help message")
          ;
  try {
    p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
    if (vars["help"].as<bool>()) {
        printHelp(desc);
        return 0;
    }
    p_o::notify(vars);

    string input_file = vars["input"].as<string>();
    string output_file = vars["output"].as<string>();
    string forcefield = vars["forcefield"].as<string>();
    string lackfile = vars["lack-file"].as<string>();
    string rtpout = vars["rtp-output"].as<string>();

    tpp::DbBase::Settings      baseSettings;
    tpp::AtomDefiner::Settings atomSettings;
    tpp::BondDefiner::Settings bondSettings;

    bool verbose = vars["verbose"].as<bool>();

    bondSettings.verbose = verbose;
    bondSettings.noqalculate = vars.count("nocalculate") == 1;
    bondSettings.ffName = forcefield;

    atomSettings.maxbonds = vars.count("max-bonds") == 1;
    atomSettings.maxdihedrals = atomSettings.maxbonds;
    atomSettings.maxangles = atomSettings.maxbonds;

    baseSettings.host     = vars["sqlserver"].as<string>();
    baseSettings.user     = vars["sqluser"].as<string>();
    baseSettings.password = vars["sqlpassword"].as<string>();
    baseSettings.port     = vars["sqlport"].as<unsigned>();
    baseSettings.dbname   = "tppforcefield"; // TODO remove hard code

    // finish analysing
    // starting work with input and output files

    if (verbose) {
      printInfo();
    } else {
      cout << format("Starting %1$s program.\n") % "TPPMKTOP";
    }
    // INPUT analysing
    tpp::InputFormat iform;
    bfs::path if_path(input_file);
    if (! bfs::is_regular_file(if_path) ) {
        tpp::Exception e("Error in input file.");
        e.add("error","Input file '"+if_path.filename().string()+"' is not a regular file!");
        throw e;
    }
    string in_ext = if_path.has_extension() ? if_path.extension().string() : "";
    in_ext = strutil::toLower(in_ext);
    if (in_ext == ".pdb")
      iform = tpp::TPP_IF_PDB;
    else if (in_ext == ".gro")
      iform = tpp::TPP_IF_GRO;
    else if (in_ext == ".g96")
      iform = tpp::TPP_IF_G96;
    else if ((in_ext == ".log") || (in_ext == ".out"))
      iform = tpp::TPP_IF_GAMSP;
    else {
       ostringstream os;
       tpp::Exception e("Couldn't determine format of input file. ");
       os <<
              "Unknown extension: \"" << in_ext <<"\" "
              "Please specify other extension.\n";
       e.add("error", os.str());
       throw e;
    }
    if (verbose) {
       cout << tpp::in_fmt_descr(iform)<<endl;
     }

  if (bondSettings.noqalculate) {
    cout << "TPPMKTOP will try to make full-determined topology!" << endl;
  }

  // Main  program body
    tpp::Topology TOP;
    tpp::StructureIO sio(false, rtpout.size() > 0); // it seems that ignore index has no meaning here
    tpp::ResidueNameGenerator rng(input_file);
    TOP.res_name = rng.getName();
    TOP.nrexcl = 3;
    sio.loadFromFile(TOP, iform, input_file.c_str());

    // initial DB queries
    tpp::DbInfo DI(baseSettings, forcefield);
    atomSettings.ffID = DI.getFFID();
    TOP.ffinclude = DI.getFFInclude().c_str();
    TOP.ffinfo = forcefield + " revision " + DI.getFFRev();
    if (verbose) {
      cout << DI.getStatistics();
    }
    // starting program body
    tpp::AtomDefiner AD(baseSettings, atomSettings, TOP);
    AD.proceed();
    AD.atomAlign();
    tpp::BondDefiner BD(baseSettings, bondSettings, TOP);
    BD.bond_align();
    tpp::save_topology(TOP, output_file.c_str(), bondSettings.noqalculate);
    tpp::save_lack(TOP, lackfile.c_str());
    if (rtpout.size() > 0) {
      tpp::save_topology_rtp(TOP, rtpout.c_str());
    }
    cout << format("Please, correct your charges according to sum: %1$8.3f.\n") % sumcharge(TOP);
  } // of global try
  catch (const tpp::SqlException &e) {
    TPPE << "TPP_SQL_EXCEPTION: " << e.what() << endl;
    return 3;
  } catch (const tpp::DbException &e) {
    TPPE << "TPP_DB_EXCEPTION: " << e.what() << endl;
    return 2;
  }
  catch (boost::program_options::error & e) {
    cerr << format("\nTPPRENUM %1% : Error in input parameters.\n\n") % VERSION;
    cerr << desc;
    return 1;
  }
  catch (const tpp::Exception &e) {
    TPPE << "  TPP_EXCEPTION: " << e.what() << endl;
    return 2;
  }
  catch(const std::exception& e)
  {
    cerr<<"  TPP crashed with std::exception:"<<e.what()<<endl;
    return 3;
  }
  catch(...)
  {
    cerr<<"  TPP crashed with unknown type of exception!"<<endl;
    return 3;
  }

  cout << "TPPMKTOP finished normally!" << endl;
  return 0;
}

void printInfo(){
  cout << format(
            "\
**********************************************************************\n\
*   Biology faculty, Department of biophysics, Erg Research Group    *\n\
*   Moscow, Lomonosov's Moscow State University                      *\n\
*   for more info, see homepage  http://erg.biophys.msu.ru/          *\n\
*                                                                    *\n\
*   Authors:       comcon1, dr.zoidberg, piton, leela, month         *\n\
*                                                                    *\n\
*   Product:       program  TPPMKTOP-%1$-6s                          *\n\
*                                                                    *\n\
*    Utilite for generating final topology from PDB file according   *\n\
* to SMARTS-patterns in database. Also you can maximize correctness  *\n\
* of topology by using special algorythmes.                          *\n\
*                                                                    *\n\
* Configured:     %2$-19s                                *\n\
**********************************************************************\n\
\n\n")
            % PACKAGE_VERSION % CONFIGURE_CDATE;
}

/*!
 * \brief Prints program info and usage to stdout.
 */
void printHelp(p_o::options_description const&_desc) {
  cout
      << format(
          "\n\
--------------------------------*****---------------------------------\n\
                                          ERG Research Group,         \n\
                                          Department of biophysics,   \n\
                                          Biology faculty, MSU, Russia\n\
\n\
           THE PART OF TOPOLOGY PREPROCESSOR PROJECT                  \n\
        ---      (comcon1, zoidberg, piton, leela, month)      ---    \n\
  TPP version: %1$-3s, compiled at %2$-8s on GCC %3$s.\n\
  BOOST version:  %4$-8s \n\
  OpenBabel version: %5$-8s \n\
  OpenBabel share data: %6$-8s \n\
\n\
                                TPPMKTOP\n\
   Utilite for checking your structure file to be suite for  next-step\n\
programs.  Also it adapts names and  position of atoms in file to make\n\
following topology file more obvious.\n\
\n\
")        % PACKAGE_VERSION % CONFIGURE_CDATE % __VERSION__ % BOOST_LIB_VERSION
          % BABEL_VERSION % BABEL_DATADIR
      << _desc << "\n\
--------------------------------*****---------------------------------\n\
"
    << endl;
}

double sumcharge(const tpp::Topology &tp) {
  double sum = 0.0;
  for (tpp::AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end();
      ++it) {
    sum += it->charge;
  }
  return sum;
}
