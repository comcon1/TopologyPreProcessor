#include "../include/db_base.hpp"
#include "global.hpp"
#include "topio.hpp"

using namespace tpp;
using namespace ut;

void tpp_exception_translator(t_exception e) {
  cerr << "ERROR FROM: " << e["procname"] << endl;
  BOOST_MESSAGE("TPP EXCEPTION FIXED!");
  e.fix_log();
}
  // define parameters for atom_definer
t_input_params par0;
t_topology TOP;

void TEST_ATALIGN() throw (t_exception);

void TEST_BONDALIGN() throw (t_exception);

int test_main( int , char* [] ) {
  PARAM_ADD(par0, "host", "localhost");
  PARAM_ADD(par0, "dbname", "tppforcefield");
  PARAM_ADD(par0, "user", "tppuser");
  PARAM_ADD(par0, "password", "estatic");
  PARAM_ADD(par0, "port", "3306");
  PARAM_ADD(par0, "ffname", "OPLS-AA");
//  PARAM_ADD(par0, "maxbonds", "on");
//  PARAM_ADD(par0, "maxangles", "on");
//  PARAM_ADD(par0, "maxdihedrals", "on");
// PARAM_ADD(par0, "noqalculate", "on");
 PARAM_ADD(par0, "qalcfile", "lack.itp");
  load_struct(TOP, TPP_IF_PDB, "hexane.pdb");
  TOP.name = "n-HEXANE";
  TOP.res_name = "HEX";

  unit_test_monitor.register_exception_translator<t_exception>( &tpp_exception_translator );
  
  ut::test_suite *test0 = BOOST_TEST_SUITE("Aligning atom types!");
  test0->add( BOOST_TEST_CASE(&TEST_ATALIGN), 0 );
  framework::run( test0 );
  ut::test_suite *test1 = BOOST_TEST_SUITE("Aligning bond types!");
  test1->add( BOOST_TEST_CASE(&TEST_BONDALIGN), 0 );
  framework::run( test1 );

  save_topology(TOP, "result.itp");

  cout << "-------------------------------------------------------------------------" << endl;
  ut::results_reporter::detailed_report( test0->p_id );
  ut::results_reporter::detailed_report( test1->p_id );
  cout << "-------------------------------------------------------------------------" << endl;

  return 0;
}

void TEST_ATALIGN() throw (t_exception) {

  atom_definer AD(par0, TOP);
  AD.proceed();
  AD.atom_align();
}
void TEST_BONDALIGN() throw (t_exception) {

  bond_definer BD(par0, TOP);
  BD.bond_align();
}
