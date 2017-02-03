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

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/format.hpp>

#define TPP_LOADFILE_TIMELIMIT 60
#define TPP_RENUMBER_TIMELIMIT 60

namespace p_o = boost::program_options;

using boost::format;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using namespace boost::numeric;

void print_help();
void print_info();

string extension(const std::string& filename){
  string::size_type ind = filename.find(".", 0);
  if (ind == string::npos) {
       return "";
  }
  return filename.substr(ind + 1);
}

int main(int argc, char * argv[]) {
  tpp::initiate_logging("tpprenum.log");
  string progname("Execution rules for TPPRENUM ");
  progname += PACKAGE_VERSION;
  p_o::options_description desc(progname);
  p_o::variables_map vars;
  desc.add_options()
      ("input,i",
          p_o::value<std::string>()->required(),
          "Input filename (any format)")
      ("output,o",
          p_o::value<std::string>()->required(),
          "Output filename (any format)")
      ("hex,x",
          p_o::value<bool>()->default_value(false)->implicit_value(false),
          "Hexadecimal numbering of main chain")
      ("verbose,v",
          p_o::value<bool>()->default_value(false)->implicit_value(false),
          "Verbose mode")
      ("ignore-index,g",  // TODO:this should be without negation!
          p_o::value<bool>()->default_value(false)->implicit_value(false),
          "Don't ignore index")
      ("rtpoutput-file,r",
          p_o::value<bool>()->default_value(false)->implicit_value(false),
          "????") // TODO: add description
      ("help,h", "Print this message")
      ;
  try {
    p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
    p_o::notify(vars);

    if (vars.count("help"))
    {
      print_help();
      return 0;
    }

    bool verbose = vars["verbose"].as<bool>();
    bool hex_flag = vars["hex"].as<bool>();
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
      print_info();
    } else {
      TPPI << format("Starting TPPRENUM-%1$s program.") % VERSION;
    }

    // INPUT analysing
    tpp::InputFormat iform;
    tpp::OutputFormat oform;
    string in_ext = extension(input_file);
    if (in_ext == "pdb")
      iform = tpp::TPP_IF_PDB;
    else if (in_ext == "gro")
      iform = tpp::TPP_IF_GRO;
    else if (in_ext == "g96")
      iform = tpp::TPP_IF_G96;
    else if ((in_ext == "log") || (in_ext == "out"))
      iform = tpp::TPP_IF_GAMSP;
    else {
      cerr << "ERROR:\n";
      cerr << "Couldn't determine format of input file. "
              "Unknown extension: \"" << in_ext <<"\" "
              "Please specify other extension.\n";
      return 1;
    }

    if (verbose) {
      TPPI<<tpp::in_fmt_descr(iform);
    }

    // OUTPUT analysing
    string out_ext = extension(output_file);
    if (out_ext == "pdb")
      oform = tpp::TPP_OF_PDB;
    else if (out_ext == "gro")
      oform = tpp::TPP_OF_GRO;
    else if (out_ext == "g96")
      oform = tpp::TPP_OF_G96;
    else {
      TPPE << "ERROR:\n"
           << "Couldn't determine format of output file."
              "Unknown extension: \"" << out_ext <<"\""
              " Please specify other extension.";
      return 1;
    }

    if (verbose) {
      TPPI<< tpp::out_fmt_descr(oform);
    }
    //
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
          topology.atoms = rnr.molRenumber(topology.atoms, tail1, hex_flag);
          } );

    io.saveToFile(topology, oform, output_file.c_str());
  } // of global try
  catch (const boost::program_options::error & e) {
    cerr << format("\nTPPRENUM %1% : Error in input parameters.\n\n") % VERSION;
    cerr << desc;
    return 1;
  }
  catch (const tpp::Exception &e) {
    TPPE << "TPP crashed with tpp::exception: "<< e.what();
    return 2;
  }
  catch(const std::exception& e)
  {
    TPPE<<"TPP crashed with std::exception: "<<e.what();
    return 3;
  }
  catch(...)
  {
    TPPE<<"TPP crashed with unknown type of exception!";
    return 3;
  }
  TPPI << "TPPRENUM finished normally!" << endl;
  return 0;
}

/// Prints brief information.
void print_info()
{
  cout << format ("\
   **********************************************************************\n\
   *   Biology faculty, Department of biophysics, Erg Research Group    *\n\
   *   Moscow, Lomonosov's Moscow State University                      *\n\
   *   for more info, see homepage  http://erg.biophys.msu.ru/          *\n\
   *                                                                    *\n\
   *   Authors:       comcon1, dr.zoidberg, piton                       *\n\
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
void print_help()
{
    cout << format("\n\
--------------------------------*****---------------------------------\n\
                                          ERG Research Group,         \n\
                                          Department of biophysics,   \n\
                                          Biology faculty, MSU, Russia\n\
\n\
           THE PART OF TOPOLOGY PREPROCESSOR PROJECT                  \n\
        ---      (comcon1, zoidberg, piton)       ---                 \n\
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
 USAGE: \n\
 tpprenum -i <input> -o <output> [-v]   \n\
      -i  the name of (I)nput-file, in PDB or GRO/G96 format.           \n\
      -o  the name of (O)utput-file, contained prepared structure.      \n\
      -v  (V)erbose mode, typing more information during the execution\n\
      -x  he(X)adecimal numbering of main (heavy-atom) chain          \n\
      -g  Do not i(G)nore index in PDB (vs sequental order).    \n\
      -h  print this message.                                         \n\
\n\
--------------------------------*****---------------------------------\n\
") % PACKAGE_VERSION % CONFIGURE_CDATE % __VERSION__ % BOOST_LIB_VERSION
         % BABEL_VERSION % BABEL_DATADIR << endl;
    throw 0;
}


