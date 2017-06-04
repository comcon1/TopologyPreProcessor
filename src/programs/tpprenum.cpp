/*! \file tpprenum.cpp
 *  \brief This file is used to launch the TPPRENUM, a program renumerates atoms from PDB.
 *
 * DB renumbering utility. The renumbering relies on the longest
 * bonded atomic sequence search. Atom names will be also modified.
 * TPPRENUM is used when atoms in  PDB file are posed in chaotic order.
 *
 */

#include "core.hpp"
#include "global.hpp"
#include "logger.hpp"
#include "exceptions.hpp"
#include "pdbutils.hpp"
#include "structio.hpp"
#include "async_call.hpp"
#include "strutil.hpp"
#include "process_options.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <sstream>

#define TPP_LOADFILE_TIMELIMIT 600
#define TPP_RENUMBER_TIMELIMIT 600

namespace p_o = boost::program_options;
namespace bfs = boost::filesystem;

using boost::format;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using namespace boost::numeric;

void printHelp(p_o::options_description const&);
void printInfo();

string extension(const std::string& filename){
  string::size_type ind = filename.find(".", 0);
  if (ind == string::npos) {
       return "";
  }
  return filename.substr(ind + 1);
}

int main(int argc, char * argv[]) {
  string progname("Execution rules for TPPRENUM ");
  progname += PACKAGE_VERSION;
  p_o::options_description desc(progname), mandatory("Mandatory settings"),
    optional("Optional settings");
  p_o::variables_map vars;
  mandatory.add_options()
      ("input,i",
          p_o::value<std::string>()->required(),
          "Input filename (any format)")
      ("output,o",
          p_o::value<std::string>()->required(),
          "Output filename (any format)")
          ;
  optional.add_options()
      ("base36,x",
          p_o::bool_switch()->default_value(false),
          "Numbering of heavy atoms in base36")
      ("verbose,v",
          p_o::bool_switch()->default_value(false),
          "Verbose mode")
      ("ignore-index,g",  // TODO:this should be without negation!
          p_o::bool_switch()->default_value(false),
          "Don't ignore index")
      ("rtpoutput-file,r",
          p_o::bool_switch()->default_value(false),
          "????") // TODO: add description
      ("help,h", p_o::bool_switch()->default_value(false),
          "Print detailed help message")
      ;
  desc.add(mandatory).add(optional);
  try {
    p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
    if (vars["help"].as<bool>()) {
        printHelp(desc);
        return 0;
    }
    p_o::notify(vars);

    bool verbose = vars["verbose"].as<bool>();
    tpp::initiate_logging("tpprenum.log", verbose);
    bool b36Flag = vars["base36"].as<bool>();
    bool ignore_index = !vars["ignore-index"].as<bool>(); // for some reason, this flag is inverted
    bool rtp_file = !vars["rtpoutput-file"].as<bool>();

    if (! ignore_index) {
      TPPI << "Non-ignoring of indexes is a DANGEROUS MODE!" << endl; // Double negation, confusing
    }

    string input_file =vars["input"].as<std::string>();
    string output_file = vars["output"].as<std::string>();

    //
    // finished parsing base cmd
    // starting work with input and output files
    //

    if (verbose) {
      printInfo();
    } else {
      TPPI << format("Starting TPPRENUM-%1$s program.") % VERSION;
    }
    tpp::InputFormat iform;
    // INPUT analysing
    tpp::processInput(input_file);
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

    if (verbose) {
      TPPI<<tpp::in_fmt_descr(iform);
    }

    // OUTPUT analysing
    tpp::OutputFormat oform;
    bfs::path of_path(output_file);
    string out_ext = of_path.has_extension() ? of_path.extension().string() : "";
    out_ext = strutil::toLower(out_ext);
    if (out_ext == ".pdb")
      oform = tpp::TPP_OF_PDB;
    else if (out_ext == ".gro")
      oform = tpp::TPP_OF_GRO;
    else if (out_ext == ".g96")
      oform = tpp::TPP_OF_G96;
    else {
       ostringstream os;
       tpp::Exception e("Couldn't determine format of output file. ");
       os <<
              "Unknown extension: \"" << out_ext <<"\" "
              "Please specify other extension.\n";
       e.add("error", os.str());
       throw e;
    }
    tpp::processOutputWithExt(output_file, out_ext.c_str());

    if (verbose) {
      TPPI<< tpp::out_fmt_descr(oform);
    }

    // Main program body
    //
    tpp::Topology topology;
    tpp::StructureIO io(ignore_index, rtp_file);

    tpp::run_with_timeout<void>(TPP_LOADFILE_TIMELIMIT,
            [&]() { io.loadFromFile(topology, iform, input_file.c_str()); }
             );


    tpp::run_with_timeout<void>(TPP_RENUMBER_TIMELIMIT,
        [&]() {
          tpp::Renumberer rnr(topology.mol, verbose);
          std::vector<unsigned> tail1 = rnr.findLongestChain();
          topology.atoms = rnr.molRenumber(topology.atoms, tail1, b36Flag);
          } );

    io.saveToFile(topology, oform, output_file.c_str());
  } // of global try
  catch (const boost::program_options::error & e) {
    cerr << format("\nTPPRENUM %1% : Error in input parameters.\n\n") % VERSION;
    cerr << desc;
    return 1;
  }
  catch (const tpp::Exception &e) {
    TPPE << "TPP crashed with tpp::exception: ";
    TPPE << e.what();
    return 2;
  }
  catch(const std::exception& e) { // do not use TPPE here!
    cerr << "TPP crashed with std::exception: \n";
    cerr << e.what() << endl;
    return 3;
  }
  catch(...) {
    cerr << "TPP crashed with unknown type of exception!";
    return 3;
  }
  TPPI << "TPPRENUM finished normally!" << endl;
  return 0;
}

/// Prints brief information.
void printInfo()
{
  cout << format ("\
   **********************************************************************\n\
   *   Biology faculty, Department of biophysics, Erg Research Group    *\n\
   *   Moscow, Lomonosov's Moscow State University                      *\n\
   *   for more info, see homepage  http://erg.biophys.msu.ru/          *\n\
   *                                                                    *\n\
   *   Authors:       comcon1, dr.zoidberg, piton, month                *\n\
   *                                                                    *\n\
   *   Product:       program  TPPRENUM-%1$-6s                          *\n\
   *                                                                    *\n\
   *    PDB renumbering utility. The renumbering relies on the longest  *\n\
   * bonded atomic sequence search. Atom names will be also modified.   *\n\
   * You need to run TPPRENUM when atoms in your PDB file are posed in  *\n\
   * chaotic order.                                                     *\n\
   *                                                                    *\n\
   *   Modified:     %2$-19s                                *\n\
   **********************************************************************\n\
   \n\n") % PACKAGE_VERSION % CONFIGURE_CDATE;
}

/*!
 * \brief Prints program info and usage to stdout.
 */
void printHelp(p_o::options_description const&_desc)
{
    cout << format("\n\
--------------------------------*****---------------------------------\n\
                                          ERG Research Group,         \n\
                                          Department of biophysics,   \n\
                                          Biology faculty, MSU, Russia\n\
\n\
           THE PART OF TOPOLOGY PREPROCESSOR PROJECT                  \n\
        ---      (comcon1, zoidberg, piton, month)       ---          \n\
  TPP version: %1$-3s, compiled at %2$-8s on GCC %3$s.\n\
  BOOST version:  %4$-8s \n\
  OpenBabel version: %5$-8s \n\
  OpenBabel share data: %6$-8s \n\
\n\
                                TPPRENUM\n\
    PDB renumbering utility. The renumbering relies on the longest  \n\
 bonded atomic sequence search. Atom names will be also modified.   \n\
 You need to run TPPRENUM when atoms in your PDB file are posed in  \n\
 chaotic order.                                                     \n\
\n\
"   ) % PACKAGE_VERSION % CONFIGURE_CDATE % __VERSION__ % BOOST_LIB_VERSION
         % BABEL_VERSION % BABEL_DATADIR
   << _desc <<
"\n\
--------------------------------*****---------------------------------\n\
"  << endl;
}


