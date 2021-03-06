#include "global.hpp"
#include "topio.hpp"
#include "pdbutils.hpp"
using namespace tpp;
using namespace ut;


void tpp_exception_translator(t_exception e) {
  cerr << "ERROR FROM: " << e["procname"] << endl;
  BOOST_MESSAGE("TPP EXCEPTION FIXED!");
}
void TEST_ITPIN();
void TEST_PDBIN();
void TEST_LONGTAIL();

int test_main( int , char* [] ) {
  unit_test_monitor.register_exception_translator<t_exception>( &tpp_exception_translator );
  
  ut::test_suite *test_pdbin = BOOST_TEST_SUITE("Structure and topology reading from file test!");
  test_pdbin->add( BOOST_TEST_CASE(&TEST_ITPIN), 0 );

//  test_pdbin->add( BOOST_TEST_CASE(&TEST_LONGTAIL), 0 );
  framework::run( test_pdbin );

  cout << "-------------------------------------------------------------------------" << endl;
  ut::results_reporter::detailed_report( test_pdbin->p_id );
  cout << "-------------------------------------------------------------------------" << endl;

  return 0;
}

void TEST_ITPIN() {
t_topology TOP;
  try {
  load_struct(TOP, TPP_IF_GAMHESS, "ets.out");  
  load_topology(TOP, "ets.itp");  
//  save_topology(TOP, "testout.itp");
  save_struct(TOP, TPP_OF_PDB, "etso.pdb");
//  save_lack(TOP, "lackout.itp");
//  ublas::matrix<double> mtx(ublas::zero_matrix<double>(27,27));
//  load_hessian(mtx, "hesse.out");
  }
  catch(t_exception e) {
    cerr << e["error"] << ": " << e["filename"] << endl;
  }
}

void TEST_PDBIN() {
t_topology TOP;
  try {
  load_struct(TOP, TPP_IF_PDB, "test.pdb");
  }
  catch(t_exception e) {
    cerr << e["error"] << ": " << e["filename"] << endl;
  }
  load_struct(TOP, TPP_IF_PDB, "hexane.pdb");
  TOP.name = "n-HEXANE";
  TOP.res_name = "HEX";
  save_struct(TOP, TPP_OF_PDB, "hexane1.pdb");
  save_struct(TOP, TPP_OF_GAMIN, "hexane1.inp");
   
}

void TEST_LONGTAIL() {
  t_topology TOP;
  load_struct(TOP, TPP_IF_PDB, "hexane.pdb");
  TOP.name = "n-HEXANE";
  TOP.res_name = "HEX";
  boost::numeric::ublas::vector<int> tail = generate_long_tail(TOP.mol);
  BOOST_CHECK(tail.size() > 1);
  t_atom_array X = mol_renum(TOP.mol, tail);
}

   