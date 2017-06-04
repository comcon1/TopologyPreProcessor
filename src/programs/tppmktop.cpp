/*! \file tppmktop.cpp
 *
 *	\brief This file provides executable for TPPMKTOP utility.
 *
 */

#include "core.hpp"
#include "global.hpp"
#include "exceptions.hpp"
#include "logger.hpp"
#include "topwriter.hpp"
#include "structio.hpp"
#include "db_base.hpp"
#include "tppnames.hpp"
#include "async_call.hpp"
#include "atom_definer.hpp"
#include "bond_definer.hpp"
#include "strutil.hpp"
#include "process_options.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <sstream>

#define TPP_LOADFILE_TIMELIMIT 600
#define TPP_MKTOP_TIMELIMIT 600
#define TPP_SQLCON_TIMELIMIT 10

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
  string progname("Execution rules for TPPMKTOP ");
  progname = progname + VERSION;
  p_o::options_description desc(progname), mandatory("Mandatory settings"),
    optional("Optional settings"), dbopts("Database settings");
  p_o::variables_map vars;

  mandatory.add_options()
        ("input,i",
            p_o::value<std::string>()->required(),
            "Input filename (any format)")
        ("output,o",
            p_o::value<std::string>()->required(),
            "Output filename (itp format)")
        ("forcefield,f",
            p_o::value<std::string>()->required(),
            "Forcefield name");
  optional.add_options()
        ("rtp-output,r",
            p_o::value<std::string>()->default_value(""),
            "Output filename (rtp format)")
        ("lack-file,l",
            p_o::value<std::string>()->default_value("lack.itp"),
            "Topology lack filename (default 'lack.itp')")
        ("expanded",
            p_o::value<bool>()->default_value(false)->implicit_value(false),
            "Create expanded topology (do not required FF includes)")
        ("separate",
            p_o::bool_switch()->default_value(false),
            "Create separated FF file <filename>_ff.itp")
        ("finalize",
            p_o::value<bool>()->default_value(false)->implicit_value(false),
            "Create final topology (don't create lack-file, overwrite dihedrals with pairs)")
        ("max-bonds,m",
            p_o::bool_switch()->default_value(false),
            "Maximize amount of bonds, angles and dihedrals by selecting other atom-types.")
        ("verbose,v", p_o::bool_switch()->default_value(false),
            "Verbose mode")
        ("help,h", p_o::bool_switch()->default_value(false),
          "Print detailed help message")
          ;
  dbopts.add_options()
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
        ("sqldb",
            p_o::value<std::string>()->default_value("tppforcefield"),
            "Mysql database (default 'tppforcefield')")
          ;
    desc.add(mandatory).add(optional).add(dbopts);
  try {
    p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
    if (vars["help"].as<bool>()) {
        printHelp(desc);
        return 0;
    }
    p_o::notify(vars);

    string input_file = vars["input"].as<string>();
    string output_file = vars["output"].as<string>();
    string output_ff("");
    string forcefield = vars["forcefield"].as<string>();
    string lackfile = vars["lack-file"].as<string>();
    string rtpout = vars["rtp-output"].as<string>();

    tpp::DbBase::Settings      baseSettings;
    tpp::AtomDefiner::Settings atomSettings;
    tpp::BondDefiner::Settings bondSettings;
    tpp::TopologyWriter::Settings twSettings;

    bool verbose = vars["verbose"].as<bool>();
    tpp::initiate_logging("tppmktop.log", verbose);

    bondSettings.verbose = verbose;
    bondSettings.expanded = vars["expanded"].as<bool>();
    bondSettings.finalize = vars["finalize"].as<bool>();

    atomSettings.verbose  = verbose;
    atomSettings.maxbonds = vars["max-bonds"].as<bool>();
    atomSettings.maxdihedrals = atomSettings.maxbonds;
    atomSettings.maxangles = atomSettings.maxbonds;

    baseSettings.host     = vars["sqlserver"].as<string>();
    baseSettings.user     = vars["sqluser"].as<string>();
    baseSettings.password = vars["sqlpassword"].as<string>();
    baseSettings.port     = vars["sqlport"].as<unsigned>();
    baseSettings.dbname   = vars["sqldb"].as<string>();

    twSettings.ffSeparate = vars["separate"].as<bool>();

    // finish analysing
    // starting work with input and output files

    if (verbose) {
      printInfo();
    } else {
      TPPI << format("Starting %s %s program.") % "TPPMKTOP" % PACKAGE_VERSION;
      cout << "------------------------------------------\n" << endl;
    }
    TPPI << "  == Verifying input & output ==";
    // INPUT analysing
    tpp::processInput(input_file);
    tpp::InputFormat iform;
    bfs::path if_path(input_file);
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
    TPPD << tpp::in_fmt_descr(iform);

    if (bondSettings.finalize)
      TPPI << "TPPMKTOP will try to make final topology!";
    if (bondSettings.expanded)
      TPPI << "TPPMKTOP will make self-consistent topology in separate file.";

    if (rtpout.size() > 0) {
      tpp::processOutputWithExt(rtpout, ".rtp");
    }
    tpp::processOutputWithExt(output_file, ".itp");
    if (twSettings.ffSeparate) {
      output_ff = output_file;
      tpp::processOutputWithExt(output_ff, ".itp", "_ff");
    }
    tpp::processOutputWithExt(lackfile, ".itp");

    // Main  program body
    tpp::Topology TOP;
    tpp::StructureIO sio(false, rtpout.size() > 0); // it seems that ignore index has no meaning here
    TOP.nrexcl = 3;

    tpp::run_with_timeout<void>(TPP_LOADFILE_TIMELIMIT,
            [&]() {
      sio.loadFromFile(TOP, iform, input_file.c_str());
    } );

    tpp::ResidueNameGenerator rng(if_path.filename().stem().string());
    TOP.res_name = rng.getName();
    TPPD << ("Using residue name: " + TOP.res_name);

    // initial DB queries
    // @todo: check connection ..
    tpp::DbInfo *DI;
    tpp::run_with_timeout<void>(TPP_SQLCON_TIMELIMIT,
        [&]() {
      DI = new tpp::DbInfo(baseSettings, forcefield);
        });
    atomSettings.ffID = DI->getFFID();
    bondSettings.ffID = DI->getFFID();
    TOP.ffinclude = DI->getFFInclude().c_str();
    TOP.ffinfo = forcefield + " revision " + DI->getFFRev();
    TOP.ffdefaults = DI->getFFDefaults();
    TPPD << ("Force field defaults: "+TOP.ffdefaults);
    TPPD << DI->getStatistics();
    delete DI;

    tpp::AtomDefiner *AD; // @TODO: change to auto_ptr
    tpp::BondDefiner *BD;

    // starting program body
    tpp::run_with_timeout<void>(TPP_MKTOP_TIMELIMIT,
        [&]() {
          AD = new tpp::AtomDefiner(baseSettings, atomSettings, TOP);
          AD->proceed();
          AD->atomAlign();
          BD = new tpp::BondDefiner(baseSettings, bondSettings, TOP);
          BD->bondAlign();
        } );
    delete AD;
    delete BD;

    // TODO: finalize & expanded
    tpp::TopologyWriter tio(twSettings);
    tio.saveITP(TOP, output_file.c_str());
    tio.saveAbsentParametersITP(TOP, lackfile.c_str());
    if (rtpout.size() > 0) {
      tio.saveRTP(TOP, rtpout.c_str());
    }
    TPPI << format("Please, correct your charges according to sum: %1$8.3f.\n") % sumcharge(TOP);

  } // of global try
  catch (const tpp::SqlException &e) {
    TPPE << "TPP_SQL_EXCEPTION: " << e.what() << endl;
    return 3;
  } catch (const tpp::DbException &e) {
    TPPE << "TPP_DB_EXCEPTION: " << e.what() << endl;
    return 2;
  }
  catch (boost::program_options::error & e) {
    cerr << format("\nTPPMKTOP %1% : Error in input parameters.\n\n") % VERSION;
    cerr << desc;
    return 1;
  }
  catch (const tpp::Exception &e) {
    TPPE << "\n  TPP_EXCEPTION: ";
    TPPE << e.what();
    return 2;
  }
  catch(const std::exception& e) { // do not use TPPE here!
    cerr << "  TPP crashed with std::exception:\n";
    cerr << e.what() << endl;
    return 3;
  }
  catch(...) {
    cerr << "  TPP crashed with unknown type of exception!" << endl;
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
\n\
Utility for automatic generation of molecular topologies. \n\
TPPMKTOP attributes atom types according to their chemical \n\
environment using database of SMARTS patterns and also assignes \n\
parameters of bonded interactions according to the force field used.\n\
\n\
TPPMKTOP takes as input molecular structure file (-i option)\n\
and the name of the force field to use (-f option). Current version \n\
of the database works with OPLS-AA force field. The input structure \n\
file should be preprocessed using TPPRENUM utility. The output \n\
topology file is written according to the .itp format (-o option).\n\
TPPMKTOP can produce .rtp topology files (-r option).\n\
\n\
Parameters of bonded inerations that were not found in the force\n\
field are summarized in the lack file (-l option).\n\
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
