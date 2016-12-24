#include "global.hpp"
#include "topio.hpp"
#include "hessian_accept.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>

namespace p_o = boost::program_options;
using tpp::cmdline;
using tpp::t_input_param;
using tpp::PARAM_ADD;
using tpp::PARAM_DEL;
using tpp::PARAM_READ;
using tpp::PARAM_EXISTS;
using boost::format;
using std::string;

void helpscreen();

int main(int argc, char * argv[]) {
 string progname =  string("Execution rules for TPPREADBONDS ") + PACKAGE_VERSION;
 p_o::options_description desc(progname);
 p_o::variables_map vars;
 desc.add_options()
	 ("input,i",p_o::value<std::string>(),"Input topology (itp format)")
	 ("output,o",p_o::value<std::string>(),"Output topology (itp format)")
	 ("hessian,e",p_o::value<std::string>(),"Input GAMESS-log (gamout format)")
	 ("lack-file,l",p_o::value<std::string>(),"Topology lack filename (default 'lack.itp')")
	 ("verbose,v","Verbose mode")
	 ("help,h", "Print this message")
 ;
 try {
    try { // u are ИДИОТ! boost::program_options
 	p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
	p_o::notify(vars);

        // boolean options
    if ( vars.count("verbose") > 1 )
      throw 1;
    PARAM_ADD(cmdline, "verbose_flag", vars.count("verbose") ? "on" : "off" );

    if (vars.count("help") == 1) helpscreen();
        // main string options
	if (vars.count("input") == 1) {
		PARAM_ADD(cmdline, "input_file", vars["input"].as<std::string>() );
        } else throw 1;
	if (vars.count("output") == 1) {
		PARAM_ADD(cmdline, "output_file", vars["output"].as<std::string>() );
        } else throw 1;
        if (vars.count("hessian") == 1) {
                PARAM_ADD(cmdline, "hessian_file", vars["hessian"].as<std::string>() );
        } else throw 1;
        if (vars.count("lack-file") == 1) {
		PARAM_ADD(cmdline, "lack_file", vars["lack-file"].as<std::string>() );
        } else throw 1;
    }
    catch (const boost::program_options::error &e) {throw 1;}
 }
 catch (int ExC) {
	 if (ExC) {
		 cerr << format("\nTPPREADBONDS %1% : Error in input parameters.\n\n") % VERSION;
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
*   Product:       program  TPPREADBONDS-%1$-6s                      *\n\
*                                                                    *\n\
*    Program reads hessian matrix from GAMESS output file and        *\n\
* calculating force constants of needed topology parameters.         *\n\
*                                                                    *\n\
*   Modified:     %2$-19s                                *\n\
**********************************************************************\n\
\n\n") % PACKAGE_VERSION % CONFIGURE_CDATE;
} else {
  cout << format("Starting TPPREADBONDS-%1$s program.\n") % PACKAGE_VERSION;
}

 // program body, using modules
 try {
   tpp::t_topology TOP;
   // setting up common topology parameters
   tpp::load_struct (TOP, tpp::TPP_IF_GAMHESS, PARAM_READ(cmdline, "hessian_file").c_str() );
   tpp::load_topology(TOP, PARAM_READ(cmdline, "input_file").c_str() );
   tpp::load_lack(TOP, PARAM_READ(cmdline, "lack_file").c_str() );
 // starting program body
   ublas::matrix<double> mtx(TOP.atoms.size()*3, TOP.atoms.size()*3);
   cout << "Loading hessian..." << flush;
   tpp::load_hessian(mtx, PARAM_READ(cmdline, "hessian_file").c_str());
   cout << "OK!\n" << flush;
   cout << "Starting converting hessian into topology..." << endl;
   tpp::accept_hessian(TOP, mtx);
   cout << "Converting hessian into topology finished!\n" << flush;
   for (tpp::t_top_map::iterator it = TOP.parameters.begin();
   it != TOP.parameters.end(); ++it) {
     tpp::t_top_coord par(*it);
     par.f = (it->type == tpp::TPP_TTYPE_BON)   ? 1 :
            ( (it->type == tpp::TPP_TTYPE_ANG)   ? 1 :
            ( (it->type == tpp::TPP_TTYPE_RBDIH) ? 3 : 1 ) );
     TOP.parameters.replace(it, par);
   }

   tpp::save_topology (TOP, PARAM_READ(cmdline, "output_file").c_str() );
 }
  catch (const tpp::t_exception &e) {
   e.fix_log();
   cerr << "TPP_EXCEPTION FROM: " << e["procname"] << endl;
   cerr << "more info see in log-file." << endl;
   return 2;
 }
 cout << "TPPMKTOP finished normally!" << endl;
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
  version %1$-3s, compiled at %2$-8s on GCC %3$s.\n\
\n\
                            TPPREADBOND\n\
\n\
--------------------------------*****---------------------------------\n\
") % PACKAGE_VERSION % CONFIGURE_CDATE % __VERSION__ << endl;
 throw 0;
}
