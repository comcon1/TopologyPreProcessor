// testing dihedral derivatives calculation))

#include <global.hpp>
#include <calc.hpp>

using tpp::t_point;

void TEST_BOND() {
  t_point p1,p2;
  p1(0) = 0.0;  p1(1) = 0.0; p1(2) = 00.;
  p2(0) = 0.5;  p2(1) = 0.5; p2(2) = 0.5;
  BOOST_CHECK( fcmp(sqrt(0.75)) == tpp::__CALC_BOND(0,p1,p2) );

  tpp::__CALC_BOND(0,p1,p1);
}

void TEST_BOND_PART() {
  t_point p1,p2;
  p1(0) = 0.0;  p1(1) = 0.0; p1(2) = 0.0;
  p2(0) = 0.5;  p2(1) = 0.5; p2(2) = 0.5;
  double d[7];
  for (int i=0;i<7;i++) {
    d[i] = tpp::__CALC_BOND(i,p1,p2);
//    cout << i << ": ==" << d[i] << endl;
  }
  BOOST_CHECK(d[1] < 0);
  BOOST_CHECK(d[2] < 0);
  BOOST_CHECK(d[3] < 0);
  BOOST_CHECK(d[4] > 0);
  BOOST_CHECK(d[5] > 0);
  BOOST_CHECK(d[6] > 0);
  BOOST_CHECK(fcmp(d[1]) == d[2]);
}

void TEST_ANG() {
  t_point p1,p2,p3;
  p1(0) = 0.0;  p1(1) = 0.0; p1(2) = 0.0;
  p2(0) = 1.0;  p2(1) = 0.0; p2(2) = 0.0;
  p3(0) = 0.0;  p3(1) = 0.0; p3(2) = 1.0;
  double d[10];
  for (int i=0;i<10;i++) 
    d[i] = tpp::__CALC_ANG(i,p3,p1,p2);
  
  BOOST_CHECK( fcmp(d[0] ) == M_PI / 2 );
  BOOST_CHECK( fcmp(tpp::__CALC_ANG(0,p2,p1,p2*8)) == 0.0 );
  BOOST_CHECK( fcmp(d[1] ) <  0 );
  BOOST_CHECK( fcmp(d[2] ) == 0 );
  BOOST_CHECK( fcmp(d[3] ) == 0 );
  BOOST_CHECK( fcmp(d[4] ) >  0 );
  BOOST_CHECK( fcmp(d[5] ) == 0 );
  BOOST_CHECK( fcmp(d[6] ) >  0 );
  BOOST_CHECK( fcmp(d[7] ) == 0 );
  BOOST_CHECK( fcmp(d[8] ) == 0 );
  BOOST_CHECK( fcmp(d[9] ) <  0 );

}

void TEST_BAD_ANG() {
  t_point p1,p2,p3;
  p1(0) = 0.0;  p1(1) = 0.0; p1(2) = 0.0;
  p2(0) = 1.0;  p2(1) = 0.0; p2(2) = 0.0;
  p3(0) = 0.0;  p3(1) = 0.0; p3(2) = 1.0;
  double d[10];
  d[1] = tpp::__CALC_ANG(2,p3,p1,p3);
  d[2] = tpp::__CALC_ANG(10, p1,p2,p3);

  BOOST_CHECK( fcmp(tpp::__CALC_ANG(0,p2,p1,p3)) == tpp::__CALC_ANG(0,p3,p1,p2) );
}

void TEST_DIH() throw (tpp::t_exception) {
  t_point p1,p2,p3,p4,p5;
  p1(0) = 0.0;  p1(1) = 0.0; p1(2) = 0.0;
  p2(0) = 1.0;  p2(1) = 0.0; p2(2) = 0.0;
  p3(0) = 1.0;  p3(1) = 0.0; p3(2) = 1.0;
  p4(0) = 1.0;  p4(1) = 1.0; p4(2) = 1.0;
  p5(0) = 1.0;  p5(1) = 0.0; p5(2) = 2.0;

  double d[14];
  for(int i=0;i<13;i++) {
    d[i] = tpp::__CALC_DIH(i,p1,p2,p3,p4);
  }

  BOOST_CHECK( fcmp(d[0]) == - M_PI / 2.0 ); 
  BOOST_CHECK( fcmp(d[0]) ==  tpp::__CALC_DIH(0,p4,p3,p2,p1) );
  BOOST_CHECK( fcmp(d[ 1]) == 0 );
  BOOST_CHECK( fcmp(d[ 2]) >  0 );
  BOOST_CHECK( fcmp(d[ 3]) == 0 );
  BOOST_CHECK( fcmp(d[ 4]) == 0 );
  BOOST_CHECK( fcmp(d[ 5]) <  0 );
  BOOST_CHECK( fcmp(d[ 6]) == 0 );
  BOOST_CHECK( fcmp(d[ 7]) >  0 );
  BOOST_CHECK( fcmp(d[ 8]) == 0 );
  BOOST_CHECK( fcmp(d[ 9]) == 0 );
  BOOST_CHECK( fcmp(d[10]) <  0 );
  BOOST_CHECK( fcmp(d[11]) == 0 );
  BOOST_CHECK( fcmp(d[12]) == 0 );
  
  // more complicated example

  p1(0) = 0.1;  p1(1) =-0.1; p1(2) = 0.2;
  p2(0) = 1.2;  p2(1) = 0.1; p2(2) = 0.3;
  p3(0) = 1.1;  p3(1) =-0.2; p3(2) = 0.9;
  p4(0) = 1.3;  p4(1) = 1.7; p4(2) = 1.2;

  for(int i=0;i<13;i++) {
    d[i] = tpp::__CALC_DIH(i,p1,p2,p3,p4);
  }
  
  // testing differential with numerical
  double delta = 1e-6, r0,r1;
  p2(0) -= delta/2;
  r0 = tpp::__CALC_DIH(0,p1,p2,p3,p4);
  p2(0) += delta;
  r1 = tpp::__CALC_DIH(0,p1,p2,p3,p4);
  BOOST_CHECK(  fabs(((r1-r0)/delta) - d[4]) < 0.01 );
  p2(0) -= delta/2;

  p3(2) -= delta/2;
  r0 = tpp::__CALC_DIH(0,p1,p2,p3,p4);
  p3(2) += delta;
  r1 = tpp::__CALC_DIH(0,p1,p2,p3,p4);
  BOOST_CHECK(  fabs(((r1-r0)/delta) - d[9]) < 0.01 );
  p3(2) -= delta/2;

  // must fail !

  p1(0) = 0.0;  p1(1) = 0.0; p1(2) = 0.0;
  p2(0) = 1.0;  p2(1) = 0.0; p2(2) = 0.0;
  p3(0) = 1.0;  p3(1) = 0.0; p3(2) = 1.0;
  p5(0) = 1.0;  p5(1) = 0.0; p5(2) = 2.0;

  tpp::__CALC_DIH(0,p1,p2,p3,p5);
}


void TEST_DIH_10000() throw (tpp::t_exception) {
  t_point p1,p2,p3,p4;

  double sum;

  p1(0) = 0.1;  p1(1) =-0.1; p1(2) = 0.2;
  p2(0) = 1.2;  p2(1) = 0.1; p2(2) = 0.3;
  p3(0) = 1.1;  p3(1) =-0.2; p3(2) = 0.9;
  p4(0) = 1.3;  p4(1) = 1.7; p4(2) = 1.2;

  for (int i=0; i<100000; i++) 
    for (int j=0; j<13; j++)
      sum += tpp::__CALC_DIH(j,p1,p2,p3,p4);

}

