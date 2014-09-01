#include "global.hpp"
#include "ublasadd.hpp"
#include "boost/numeric/bindings/lapack/gesv.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"

using namespace tpp;
namespace lapack = boost::numeric::bindings::lapack;

void helpscreen();

void tpp_exception_translator(t_exception ex) {
   cerr << "EXECPTION FROM: " << ex["procname"] << endl;
   BOOST_MESSAGE("TPP EXCEPTION INSTANCE CAUGHT!");
}

void TEST_LINAL() {
  ublas::matrix<double> A(4,4);
  A(0,0) = 4; A(0,1) = 3; A(0,2) = 1; A(0,3) = 1;
  A(1,0) =-3; A(1,1) =-2; A(1,2) = 1; A(1,3) = 2;
  A(2,0) =-2; A(2,1) = 1; A(2,2) =-4; A(2,3) =-3;
  A(3,0) = 1; A(3,1) = 0; A(3,2) = 2; A(3,3) = 4;
  cout << svd_rank(A) << endl;
  ublas::matrix<double,ublas::column_major> B(3,3);
  B = ublas::subrange(A, 0, 3, 0, 3);
  ublas::matrix<double, ublas::column_major> C(3,1);
  C = ublas::subrange(A, 0, 3, 3, 4);
  cout << C << endl;
  lapack::gesv(B,C);
  cout << C << endl;
  svd_inv(A,C);
  cout << prod(C,A) << endl;
}


int test_main(int argc, char * argv[]) {

  ut::unit_test_monitor.register_exception_translator<t_exception>( &tpp_exception_translator );
  ut::test_suite *test_1 = BOOST_TEST_SUITE("LINAL CALCULUS");
  test_1->add( BOOST_TEST_CASE(&TEST_LINAL), 0);
  ut::framework::run( test_1 );

  return 0;
}

