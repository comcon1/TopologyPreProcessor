// functors with standart functions

//-------------------------------------------
// of dihedral potential energy
namespace tpp {
class symdih_functor {
  public:
  static const char c_fun   = 1;
  static const char c_param = 3;
  double param[3];
  inline double operator() (double x) {
    return param[0]*(1+cos(param[1]*x+param[2]));
  };
};

class rbdih_functor {
  public:
  static const char c_fun   = 3;
  static const char c_param = 5;
  double param[5];
  inline double operator() (double x) {
    return param[0]+param[1]*cos(x)+param[2]*cos(x)*cos(x)+
      param[3]*cos(x)*cos(x)*cos(x)+param[4]*cos(x)*cos(x)*cos(x)*cos(x)+
      param[5]*cos(x)*cos(x)*cos(x)*cos(x)*cos(x);
  };
};

class impdih_functor {
  public:
    static const char c_fun   = 2;
    static const char c_param = 2;
    double param[2];
    inline double operator(double x) {
      return 0.5*param[2]*(x-param[1])*(x-param[1]);
    };  
};

//---------------------------------
// bond functors

class g96bon_functor {
  public:
    static const char c_fun   = 2;
    static const char c_param = 2;
    double param[2];
    inline double operator(double x) {
      return 0.5*param[2]*(x-param[1])*(x-param[1])*(x-param[1])*(x-param[1]);
    };  
};

class harbon_functor {
  public:
    static const char c_fun   = 1;
    static const char c_param = 2;
    double param[2];
    inline double operator(double x) {
      return 0.5*param[2]*(x-param[1])*(x-param[1]);
    };  
};

//---------------------------------
// angle functors

class g96ang_functor {
  public:
    static const char c_fun   = 2;
    static const char c_param = 2;
    double param[2];
    inline double operator(double x) {
      return param[2]*cos(x-param[1])*cos(x-param[1]);
    };  
};

class harang_functor {
  public:
    static const char c_fun   = 1;
    static const char c_param = 2;
    double param[2];
    inline double operator(double x) {
      return 0.5*param[2]*(x-param[1])*(x-param[1]);
    };  
};

}

