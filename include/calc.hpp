 /* This C++ file created with the help of GiNaC package
  * by comcon1.
  * Erg Research Group. Apr, 2007.
  */

 // DIHEDRAL ESTIMATION..OK!
 // VALENCE ANGLE ESTIMATIONS...OK!
#ifndef TPP_CALC_H
#define TPP_CALC_H

#include "global.hpp"

namespace tpp {

extern double __CALC_DIH(int i, const Point& v1, const Point& v2, 
    const Point& v3, const Point& v4) throw (t_exception);

extern double __CALC_ANG(int i, const Point& v1,const Point& v2,const Point& v3);

extern double __CALC_BOND(int i, const Point& v1,const Point& v2);


}

#endif

