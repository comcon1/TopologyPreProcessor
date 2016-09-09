#include "topio.hpp"
#include "db_scanner.hpp"
#include "testcase.hpp"

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
double sumcharge(const tpp::t_topology &);

int main(int argc, char * argv[]) {
  string progname("Execution rules for TPPMKTEST ");
  progname = progname + VERSION;
  p_o::options_description desc(progname);
  p_o::variables_map vars;
  desc.add_options()
          ("input,i",p_o::value<std::string>(),"Input filename (PDB format)")
          ("testid,t",p_o::value<unsigned>(),"Test molecule ID (zero for all)")
          ("forcefield,f",p_o::value<std::string>(),"Forcefield name")
 
          ("sqlserver,s",p_o::value<std::string>(),"Mysql-server adress (default 'localhost')")
          ("sqlport,t",p_o::value<unsigned>(),"Mysql-server port (default '3306')")
          ("sqluser,u",p_o::value<std::string>(),"Mysql-user (default 'tppuser')")
          ("sqlpassword,p",p_o::value<std::string>(),"Mysql-password (default 'estatic')")

          ("verbose,v","Verbose mode")
          ("help,h", "Print this message")
  ;
  try {
     try { 
        // parsing boost::program_options
  	p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
        p_o::notify(vars);
     
         // boolean option
        if (vars.count("verbose") > 1)  
          throw 1;
        PARAM_ADD(cmdline, "verbose_flag", vars.count("verbose") ? "on" : "off" );
        // tests are performed in nocalculate only
        PARAM_ADD(cmdline, "nocalculate_flag",  "on" );

        if (vars.count("help") == 1) helpscreen();
        
        // main string options
        if (vars.count("input") == 1) {
        	PARAM_ADD(cmdline, "input_file", vars["input"].as<std::string>() );
        } else if (vars.count("testid") == 1) {
                unsigned o = vars["testid"].as<unsigned>();
        	PARAM_ADD(cmdline, "testid", boost::lexical_cast<string>(vars["testid"].as<unsigned>()) );
        } else throw 1;
        if (vars.count("forcefield") == 1) {
        	PARAM_ADD(cmdline, "forcefield", vars["forcefield"].as<std::string>() );
        } else throw 1;

         
        // SQL parameters
        if (vars.count("sqlserver") == 1) {
        	PARAM_ADD(cmdline, "sqlserver", vars["sqlserver"].as<std::string>() );
        } else if (vars.count("sqlserver") == 0 ) {
        	PARAM_ADD(cmdline, "sqlserver", "localhost" );
        } else throw 1;
        if (vars.count("sqluser") == 1) {
        	PARAM_ADD(cmdline, "sqluser", vars["sqluser"].as<std::string>() );
        } else if (vars.count("sqluser") == 0 ) {
        	PARAM_ADD(cmdline, "sqluser", "tppuser" );
        } else throw 1;
        if (vars.count("sqlport") == 1) {
                unsigned o = vars["sqlport"].as<unsigned>();
        	PARAM_ADD(cmdline, "sqlport", boost::lexical_cast<string>(vars["sqlport"].as<unsigned>()) );
        } else if (vars.count("sqlport") == 0 ) {
        	PARAM_ADD(cmdline, "sqlport", "3306" );
        } else throw 1;
        if (vars.count("sqlpassword") == 1) {
        	PARAM_ADD(cmdline, "sqlpassword", vars["sqlpassword"].as<std::string>() );
        } else if (vars.count("sqlpassword") == 0 ) {
        	PARAM_ADD(cmdline, "sqlpassword", "estatic" );
        } else throw 1;
     }
     catch (boost::program_options::error &e) {throw 1;}
  }
  catch (int ExC) {
          if (ExC) {
         	 cerr << format("\nTPPMKTEST %1% : Error in input parameters.\n\n") % VERSION;
          	 cerr << desc;
          }
          return(ExC);
  }
  // finish analysing
  // starting work with input and output files
  
  if (PARAM_READ(cmdline, "verbose_flag") == "on") {
  cout << format ("\
**********************************************************************\n\
*   Biology faculty, Department of biophysics, Erg Research Group    *\n\
*   Moscow, Lomonosov's Moscow State University                      *\n\
*   for more info, see homepage  http://erg.biophys.msu.ru/          *\n\
*                                                                    *\n\
*   Authors:       comcon1, dr.zoidberg, piton, leela                *\n\
*                                                                    *\n\
*   Product:       program  TPPMKTEST-%1$-6s                         *\n\
*                                                                    *\n\
*    Utilite for generating test case from PDB file. This run produ- *\n\
* ces the same topology as `tppmktop -n` but do not show it to you   *\n\
* but write it to db or compare with already written records.        *\n\
*                                                                    *\n\
* Configured:     %2$-19s                                *\n\
**********************************************************************\n\
\n\n") % PACKAGE_VERSION % CONFIGURE_CDATE;
  } else {
    cout << format("Starting %1$s program.\n") % "TPPMKTEST";
  }

  // INPUT analysing
  tpp::t_iformat iform;
  string::size_type ind = PARAM_READ(cmdline, "input_file").find(".",0);
  if ( ind == string::npos) {
    cerr << "ERROR:\n";
    cerr << "Couldn't determine format of input file. Please specify extension.\n";
    return 1;
  }
  string subs = PARAM_READ(cmdline, "input_file").substr(ind+1);
  if (subs == "pdb") iform = tpp::TPP_IF_PDB;
  else if ( (subs == "gro") || (subs == "g96") || (subs == "out") ) {
    cerr << "ERROR:\n";
    cerr << "Only PDB format should be used for writing test case.\n";
    return 1;
  } else {
    cerr << "ERROR:\n";
    cerr << "Couldn't determine format of input file. Please specify other extension.\n";
    return 1;
  }

  if (PARAM_READ(cmdline, "verbose_flag") == "on") {
      cout << "Input file format: Protein Data Bank." << endl; 
  }

  // program body, using modules
  try {
    tpp::t_topology TOP;

    // setting up common topology parameters
    TOP.res_name = PARAM_READ(cmdline, "input_file").substr(0,3);
    TOP.nrexcl = 3;
    // ;-)
    tpp::load_struct_fname (TOP, iform, PARAM_READ(cmdline, "input_file").c_str() );
    // customization of 2-nd level parameters
    tpp::t_input_params par0;
    PARAM_ADD(par0, "host", PARAM_READ(cmdline, "sqlserver") );
    PARAM_ADD(par0, "dbname", "tppforcefield" );
    PARAM_ADD(par0, "user", PARAM_READ(cmdline, "sqluser"));
    PARAM_ADD(par0, "password", PARAM_READ(cmdline, "sqlpassword") );
    PARAM_ADD(par0, "port", PARAM_READ(cmdline, "sqlport"));
    PARAM_ADD(par0, "ffname", PARAM_READ(cmdline, "forcefield") );
    // initial DB queries
    tpp::db_info DI(par0);
    PARAM_ADD(par0, "ffid", boost::lexical_cast<string>(DI.get_ffid()) );
    TOP.ffinclude = DI.get_ffinclude().c_str();
    TOP.ffinfo = PARAM_READ(par0, "ffname") + " revision " + DI.get_ffrev();
    
    // verbose stat output

    if (PARAM_READ(cmdline, "verbose_flag") == "on") {
      cout << DI.get_statistics();
    }

    // preparing the topology

    tpp::atom_definer AD(par0, TOP);
    AD.proceed();
    AD.atom_align();
    tpp::bond_definer BD(par0, TOP);
    BD.bond_align();

    // performs checks over the topology

    double charge = sumcharge(TOP);
    
    cout << format("Processing molecule with total charge: %1$8.3f.\n") % charge;

    if ( abs(charge - (int) charge) > 1e-5) {
      tpp::t_input_params params;
      PARAM_ADD(params, "procname", "tppmktest::main");
      PARAM_ADD(params, "error", "Error in charge check.");
      throw tpp::t_exception("Your charge is not integer!", params);
    }

    // -- WRITE TESTCASE mode --
    cout << endl;
    
    if (PARAM_EXISTS(cmdline, "input_file")) {
      cout << format("Generated topology will be written as a test case.\n") << endl;

      tpp::db_testcase dtc(par0);
      int molid = dtc.write_testmolecule(PARAM_READ(cmdline, "input_file").c_str(), (int)charge, "Generated by TPPMKTEST");

      cout << format("[TC] Molecule was inserted with ID: %1$d.") % molid << endl;

      dtc.write_testcase(TOP, molid);

    }

    // -- PERFORM TEST mode --

  }
  catch (tpp::t_sql_exception &e) {
    e.fix_log();
    cerr << "TPP_SQL_EXCEPTION FROM: " << e["procname"] << endl;
    cerr << "more info see in log-file." << endl;
    return 3;
  }
  catch (tpp::t_db_exception &e) {
    e.fix_log();
    cerr << "TPP_DB_EXCEPTION FROM: " << e["procname"] << endl;
    cerr << "more info see in log-file." << endl;
    return 2;
  }
  catch (tpp::t_exception &e) {
    e.fix_log();
    cerr << "TPP_EXCEPTION FROM: " << e["procname"] << endl;
    cerr << "more info see in log-file." << endl;
    return 1;
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
        ---      (comcon1, zoidberg, piton, leela)      ---           \n\
  TPP version: %1$-3s, compiled at %2$-8s on GCC %3$s.\n\
  BOOST version:  %4$-8s \n\
  OpenBabel version: %5$-8s \n\
  OpenBabel share data: %6$-8s \n\
\n\
                                TPPMKTEST\n\
                                "
                  //TODO: write correct help message!
                                "\
\n\
 USAGE: \n\
 tppmktest [ -i <inp> | -t <id> ] -f <f.field> [other opt-s]\n\
\n\
      -i  the name of (I)nput-file, in PDB or GRO/G96 format.           \n\
      -t  test the ID of molecule in `selftest` table                   \n\
            (use zero for test all molecules from DB).                  \n\
      -f  the (F)orcefield name (f.i. OPLS-AA)                          \n\
      -v  (V)erbose mode, typing more information during the execution\n\
 [ database options ] \n\
      -s  MySQL (S)erver host name or IP\n\
      -t  MySQL server (P)ort number\n\
      -u  MySQL (U)ser\n\
      -p  MySQL (P)assword\n\
      -h  print this message.                                         \n\
\n\
--------------------------------*****---------------------------------\n\
") % PACKAGE_VERSION % CONFIGURE_CDATE % __VERSION__ % BOOST_LIB_VERSION 
   % BABEL_VERSION % BABEL_DATADIR << endl;
 throw 0;
}

double sumcharge(const tpp::t_topology &tp)  {
      double sum = 0.0;
      for (tpp::t_atom_array::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
            sum += it->charge;
      }
      return sum;
}
