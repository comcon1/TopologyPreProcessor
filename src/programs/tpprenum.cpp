#include "global.hpp"
#include "pdbutils.hpp"
#include "topio.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>

namespace p_o = boost::program_options;
using tpp::cmdline;
using tpp::PARAM_ADD;
using tpp::PARAM_DEL;
using tpp::PARAM_READ;
using tpp::PARAM_EXISTS;
using tpp::t_input_param;
using boost::format;
using std::string;
void helpscreen();

int main(int argc, char * argv[]) {
 string progname("Execution rules for TPPRENUM "); progname += VERSION;
 p_o::options_description desc(progname);
 p_o::variables_map vars;
 desc.add_options()
	 ("input,i",p_o::value<std::string>(),"Input filename (any format)")
	 ("output,o",p_o::value<std::string>(),"Output filename (any format)")
         ("hex,x","Hexadecimal numbering of main chain")
	 ("verbose,v","Verbose mode")
	 ("ignore-index,g","Don't ignore index")
	 ("help,h", "Print this message")
 ;
 try {
    try { // блок вылавливания исключений boost::program_options
 	p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
	p_o::notify(vars);
    if ( (vars.count("verbose") > 1) || (vars.count("ignore-index") > 1) ) throw 1;

    PARAM_ADD(cmdline, "verbose_flag", vars.count("verbose") ? "on" : "off" );
    PARAM_ADD(cmdline, "hex_flag",     vars.count("hex") ? "on" : "off" );
    PARAM_ADD(cmdline, "ignore_index", vars.count("ignore-index") ? "off" : "on" );
    
        if (vars.count("ignore-index") == 1) {
            cout << "Non-ignoring of indexes is a DANGEROUS MODE!" << endl;
        }
	if (vars.count("help") == 1) helpscreen();
	if (vars.count("input") == 1) {
		PARAM_ADD(cmdline, "input_file", vars["input"].as<std::string>() );
        }
	else throw 1;
	if (vars.count("output") == 1) {
		PARAM_ADD(cmdline, "output_file", vars["output"].as<std::string>() );
        }
	else throw 1;
    }
    catch (boost::program_options::error &e) {throw 1;}
 }
 catch (int ExC) {
	 if (ExC) {
		 cerr << format("\nTPPRENUM %1% : Error in input parameters.\n\n") % VERSION;
	 	 cerr << desc;
	 }
         return(ExC);
 }
// finish analysing
// starting work with input and output files
//
 
if (PARAM_READ(cmdline, "verbose_flag") == "on") {
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
*    PDB renumbering utility. The renumbering relise on the longest  *\n\
* bonded atomic sequence search. Atom names will be also modified.   *\n\
* You need to run TPPRENUM when atoms in your PDB file are posed in  *\n\
* chaotic order.                                                     *\n\
*                                                                    *\n\
*   Modified:     %2$-19s                                *\n\
**********************************************************************\n\
\n\n") % VERSION % BUILD_DATE;
} else {
  cout << format("Starting TPPRENUM-%1$s program.\n") % VERSION;
}



 // INPUT analysing
 tpp::t_iformat iform;
 tpp::t_oformat oform;
 string::size_type ind = PARAM_READ(cmdline, "input_file").find(".",0);
 if ( ind == string::npos) {
   cerr << "ERROR:\n";
   cerr << "Couldn't determine format of input file. Please specify extension.\n";
   return 1;
 }
 string subs = PARAM_READ(cmdline, "input_file").substr(ind+1);
      if (subs == "pdb") iform = tpp::TPP_IF_PDB;
 else if (subs == "gro") iform = tpp::TPP_IF_GRO;
 else if (subs == "g96") iform = tpp::TPP_IF_G96;
 else if ( (subs == "log") || (subs == "out") ) iform = tpp::TPP_IF_GAMSP;
 else {
   cerr << "ERROR:\n";
   cerr << "Couldn't determine format of input file. Please specify other extension.\n";
   return 1;
 }

 if (PARAM_READ(cmdline, "verbose_flag") == "on") {
  switch (iform) {
    case   tpp::TPP_IF_PDB: cout << "Input file format: Protein Data Bank." << endl; break;
    case   tpp::TPP_IF_GRO: cout << "Input file format: GROmacs structure." << endl; break;
    case   tpp::TPP_IF_G96: cout << "Input file format: Gromos 96 structure." << endl; break;
    case tpp::TPP_IF_GAMSP: cout << "Input file format: GAMess output." << endl; break;
  };
 }

 // OUTPUT analysing
 ind = PARAM_READ(cmdline, "output_file").find(".",0);
 if ( ind == string::npos) {
   cerr << "ERROR:\n";
   cerr << "Couldn't determine format of output file. Please specify extension.\n";
   return 1;
 }
 subs = PARAM_READ(cmdline, "output_file").substr(ind+1);
      if (subs == "pdb") oform = tpp::TPP_OF_PDB;
 else if (subs == "gro") oform = tpp::TPP_OF_GRO;
 else if (subs == "g96") oform = tpp::TPP_OF_G96;
 else {
   cerr << "ERROR:\n";
   cerr << "Couldn't determine format of output file. Please specify other extension.\n";
   return 1;
 }

 if (PARAM_READ(cmdline, "verbose_flag") == "on") {
  switch (oform) {
    case   tpp::TPP_OF_PDB: cout << "Output file format: Protein Data Bank." << endl; break;
    case   tpp::TPP_OF_GRO: cout << "Output file format: GROmacs structure." << endl; break;
    case   tpp::TPP_OF_G96: cout << "Output file format: Gromos 96 structure." << endl; break;
  };
 }

 // program body, using modules
 try{
   tpp::t_topology TOP;
   tpp::load_struct (TOP, iform, PARAM_READ(cmdline, "input_file").c_str() );
   ublas::vector<unsigned> tail1 = tpp::generate_long_tail1(TOP.mol);
   TOP.atoms = tpp::mol_renum1(TOP.mol, TOP.atoms, tail1 );
   tpp::save_struct (TOP, oform, PARAM_READ(cmdline, "output_file").c_str() ); 
 } catch (tpp::t_exception e) {
   cerr << "  TPP_EXCEPTION FROM: " << e["procname"] << endl;
   cerr << "  With following error: " << e["error"] << endl;
   cerr << "  more info see in log-file." << endl;
   return 2;
 }
 cout << "TPPRENUM finished normally!" << endl;
 return 0;
}

void helpscreen()
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
    PDB renumbering utility. The renumbering relise on the longest  \n\
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
") % PACKAGE_VERSION % CONFIGURE_CDATE % _VERSION % BOOST_LIB_VERSION 
   % BABEL_VERSION % BABEL_DATADIR << endl;
 throw 0;
}


