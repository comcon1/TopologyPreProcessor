#include <cmath>

class testclass {
 public:
  static const char c_param = 3;
  double param[3];
  inline double operator() (double x) {
   return param[0]*(1+cos(param[1]*x+param[2]));
 }
};

class testclass1 {
 public:
  static const char c_param = 3;
  double param[3];
  inline double operator() (double x) {
   return param[0]+(1+cos(param[1]*x+param[2]));
 }
};


      