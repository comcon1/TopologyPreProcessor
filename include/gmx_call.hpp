#ifndef TPP_GMXCALL_H
#define TPP_GMXCALL_H


//       - GROMACS -
#ifdef TPP_GMX_EXTERNAL
 #define TPP_MDRUN_EXEC "mdrun"
 #define TPP_GROMPP_EXEC "grompp"
 #define TPP_G_NMEIG_EXEC "g_nmeig"
 #define TPP_GMXLIB "/usr/share/gromacs/top"
#endif

#include "gromacs/types/inputrec.h"
t_inputrec standrec;

#endif
 