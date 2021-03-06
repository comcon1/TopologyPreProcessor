#include "global.hpp"
#include "calc.hpp"
#include "test_calc.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>

namespace p_o = boost::program_options;
using namespace tpp;

void helpscreen();

void tpp_exception_translator(t_exception ex) {
   cerr << "EXECPTION FROM: " << ex["procname"] << endl;
   BOOST_MESSAGE("TPP EXCEPTION INSTANCE CAUGHT!");
}


int test_main(int argc, char * argv[]) {

 char *progname = (char*)malloc(sizeof(char)*400);
 snprintf(progname, 400, "Execution rules for %s", PROG_NAME);
 p_o::options_description desc(progname);
 p_o::variables_map vars;
 desc.add_options()
	 ("param,p",p_o::value<string>(),"Sample string argument")
	 ("number,n",p_o::value<unsigned>(),"Sample numberic-type argument")
	 ("verbose,v","Verbose mode")
	 ("help,h", "Print this message")
 ;
 try {
    try { 
        p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
	p_o::notify(vars);
    if (vars.count("verbose") > 1) throw 1;
    PARAM_ADD(cmdline, "verbose", vars.count("verbose") ? "true" : "false");  
   if (vars.count("help") == 1) {
     helpscreen();
     return 0;
   }
   if (vars.count("number") == 1)
     PARAM_ADD(cmdline, "number", lexical_cast<string>(vars["number"].as<unsigned>()).c_str() );
   else throw 1;
   if (vars.count("param") == 1) 
     PARAM_ADD(cmdline, "param", vars["param"].as<string>().c_str());
   else throw 1;
    }
    catch (boost::program_options::error &e) {throw 1;}
 }
 catch (int ExC) {
	 if (ExC) {
		 cerr << format("\n%1% : Error in input parameters.\n\n") % PROG_NAME;
	 	 cerr << desc;
	 }
         return(ExC);
 }

  ut::unit_test_monitor.register_exception_translator<t_exception>( &tpp_exception_translator );
  ut::test_suite *test_bond = BOOST_TEST_SUITE("BOND CALCULATIONS");
  test_bond->add( BOOST_TEST_CASE(&TEST_BOND), 1);
  test_bond->add( BOOST_TEST_CASE(&TEST_BOND_PART), 0);
  ut::framework::run( test_bond );
  ut::test_suite *test_ang = BOOST_TEST_SUITE("ANGLE CALCULATIONS");
  test_ang->add( BOOST_TEST_CASE(&TEST_ANG), 0);
  test_ang->add( BOOST_TEST_CASE(&TEST_BAD_ANG), 2);
  ut::framework::run( test_ang );
  ut::test_suite *test_dih = BOOST_TEST_SUITE("DIHEDRAL CALCULATIONS");
  test_dih->add( BOOST_TEST_CASE(&TEST_DIH), 1);
  test_dih->add( BOOST_TEST_CASE( &TEST_DIH_10000 ), 1 );
  
  ptime t0 = microsec_clock::local_time();
  ut::framework::run( test_dih );
  time_duration diff = microsec_clock::local_time() - t0;

  cout << "-------------------------------------------------------------------------" << endl;
  ut::results_reporter::detailed_report( test_bond->p_id );
  ut::results_reporter::detailed_report( test_ang->p_id );
  ut::results_reporter::detailed_report( test_dih->p_id );
  cout << "Wasted " << diff << " for 100'000 full dihedrals!" << endl;
  cout << "-------------------------------------------------------------------------" << endl;

  cout << format("Your parameters: number (%1%), param(%2%)") 
    % PARAM_READ(cmdline,"number").c_str() % PARAM_READ(cmdline,"param").c_str() << endl;

  runtime.log_write("\n---\nProgram body logged here\n---\n");
  
  cout << "TRYING TO SERIALIZE MAP" << endl;

  return 0;
}

void helpscreen() {
  cout << "\
    ******************************************************\n\
\n\
              HELP SCREEN FOR TPP PROGRAM\n\
\n\
    ******************************************************" << endl;
}
       