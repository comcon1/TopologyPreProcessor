#include "test.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;
using namespace std;

struct struct_t {
  double a,b;
  int i,j;
  struct_t(double a0, double b0, int i0, int j0): 
  a(a0),b(b0),i(i0),j(j0){;};
  struct_t(): a(0), b(0), i(0), j(0) {;};
};

int main() {
  testclass  c;
  testclass1 *c1, *c2;
  char serial[sizeof(testclass1)];
  c1 = new testclass1();
  c.param[0] = 1.0;
  c.param[1] = -1.0;
  c.param[2] = -2;
  c1->param[0] = 120.1;
  c1->param[1] = -01.0;
  c1->param[2] = -2;

  cout << "start - "  << (*c1)(1.5) << endl;
  cout << "serialization " << sizeof(testclass1) << endl;
  char* i0 = (char *)c1;
  char* i = i0;
  while ( (long long)i < (long long)i0 + sizeof(testclass1) ) {
//    cout << i;
    printf("%8d : %4d\n ", i-i0 ,(short)(*i));
    i++;
    serial[i-i0] = *i;
  }
  cout << "finished" << endl; 
  delete c1;
  c2 = new testclass1(); 
  i0 = (char*)c2;
  i = i0;
  cout << "before recover - " << (*c2)(1.5) << endl;
  while ( (long long)i < (long long)i0 + sizeof(testclass1) ) {
    *i = serial[i-i0];
    i++;
  }
  cout << "test recover - " << (*c2)(1.5) << endl;
  ublas::vector<double> vec1(3), vec2(3);
  vec1(0) = 1.0;
  vec1(1) = -0.1;
  vec1(2) = 0.2342;
  vec2(0) = 6.0;
  vec2(1) = -3.1;
  vec2(2) = 2.2342;
  cout << vec1 << endl;
  cout << vec2 << endl;
  cout <<  inner_prod(vec1, vec2) << endl;
  cout << vec1(0)*vec2(0) + vec1(1)*vec2(1) + vec1(2)*vec2(2) << endl;
  cout << norm_1(vec1) << endl;
  cout << norm_2(vec1) << endl;
  cout << vec1(0)*vec1(0) + vec1(1)*vec1(1) + vec1(2)*vec1(2) << endl;
  cout << sqrt(vec1(0)*vec1(0) + vec1(1)*vec1(1) + vec1(2)*vec1(2)) << endl;
  struct_t st;
  st.a = 345.34;
  double *dptr = &(st.a);
  cout << *dptr << endl;
  *dptr = 200;
  cout << st.a << endl;
  vector<struct_t> pe;
  pe.push_back(struct_t(1.01,2.02,3,4));
  pe.push_back(struct_t(-1.01,-2.02,-3,-4));
  dptr = &(pe[0].a);
  cout << pe[0].a << endl;
  *dptr = 152.0;
  cout << pe[0].a << endl;

  return 0;
}
 