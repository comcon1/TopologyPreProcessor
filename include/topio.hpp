#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "global.hpp"

namespace tpp {

static const char * top_comment = 
"; TPP - topology generator \n"
"; created by Erg Research Group\n"
"; MSU, Biology Faculty, Department of Biophysics\n"
"; ----------------------------------------------\n";

extern void save_topology(t_topology &, const char *) throw (t_exception);

extern void save_topology_rtp(t_topology &, const char *) throw (t_exception);

extern void load_topology(t_topology &, const char *) throw (t_exception);

extern void load_lack(t_topology &, const char *) throw (t_exception);

extern void check_topology(t_topology &) throw (t_exception);

extern void load_struct(t_topology &, t_iformat, const char *) throw (t_exception);

extern void save_struct(t_topology &, t_oformat, const char *) throw (t_exception);

extern void save_lack(t_topology &, const char *);

#if ENABLE_GAMESS_FEATURES
extern void load_hessian(ublas::matrix<double>&, const char *) throw (t_exception);
#endif

//TODO: also serialization should be included

}
#endif

