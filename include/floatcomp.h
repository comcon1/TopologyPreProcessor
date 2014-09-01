#ifndef _FLOAT_COMPARE
#define _FLOAT_COMPARE

const double flt_epsilon = 0.00000001;

#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric;

template <class _Val,class _Prec>
class dcmp
{
public:
   typedef _Val  value_type;
   typedef _Prec prec_type;

public:
   template <class _Tp>
   dcmp(const _Tp& val): m_Val(val) {}
   bool operator== (const dcmp& rhs) const { return fabs (m_Val - rhs.m_Val) < prec_type::eps(); }
   bool operator!= (const dcmp& rhs) const { return !(*this == rhs) ; }
   bool operator<  (const dcmp& rhs) const { return *this != rhs && m_Val < rhs.m_Val ; }
   bool operator>  (const dcmp& rhs) const { return *this != rhs && m_Val > rhs.m_Val ; }
   bool operator<= (const dcmp& rhs) const { return *this == rhs || m_Val < rhs.m_Val ; }
   bool operator>= (const dcmp& rhs) const { return *this == rhs || m_Val > rhs.m_Val ; }
   friend bool operator== (const value_type& Val, const dcmp& rhs) { return dcmp(Val) == rhs ; }
   friend bool operator!= (const value_type& Val, const dcmp& rhs) { return dcmp(Val) != rhs ; }
   friend bool operator<  (const value_type& Val, const dcmp& rhs) { return dcmp(Val) < rhs ;  }
   friend bool operator>  (const value_type& Val, const dcmp& rhs) { return dcmp(Val) > rhs ;  }
   friend bool operator<= (const value_type& Val, const dcmp& rhs) { return dcmp(Val) <= rhs ; }
   friend bool operator>= (const value_type& Val, const dcmp& rhs) { return dcmp(Val) >= rhs ; }

private:
   value_type m_Val;
};

template <class _Val, int _n, class _Prec>
class dcmp_v {
  public:
    typedef ublas::bounded_vector<_Val,_n> value_type;
    typedef _Prec prec_type;
    template <class _Tp>
      dcmp_v(const _Tp& val): m_Val(val) {}
    dcmp_v(const value_type& val): m_Val(val) {}
    bool operator== (const dcmp_v& rhs) const { 
      for (int i=0;i<_n;i++)
        if  ( fabs (m_Val(i) - rhs.m_Val(i)) > prec_type::eps() )
          return false;
      return true;
    }
    bool operator!= (const dcmp_v& rhs) const { return !(*this == rhs) ; }
    friend bool operator== (const value_type& Val, const dcmp_v& rhs) { return dcmp_v(Val) == rhs ; }
    friend bool operator!= (const value_type& Val, const dcmp_v& rhs) { return dcmp_v(Val) != rhs ; }
  private:
    value_type m_Val;
};

struct _FloatPrec {
    static double eps() { return flt_epsilon; }
};

// for bounded_vector types
template<class _Val, std::size_t _n>
   inline dcmp_v<_Val, _n, _FloatPrec> vcmp(const ublas::bounded_vector<_Val, _n>& val) { 
     return dcmp_v<_Val, _n, _FloatPrec>(val); 
   }

// for numeric types
template<class _Val>
   inline dcmp<_Val,_FloatPrec> fcmp(const _Val& val) { 
     return dcmp<_Val,_FloatPrec>(val); 
   }

#endif

