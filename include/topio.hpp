#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "global.hpp"

namespace tpp {

static const char * top_comment =
"; ----------------------------------------------\n" 
"; TPP - topology generator version " PACKAGE_VERSION " \n"
"; created by Erg Research Group\n"
"; MSU, Biology Faculty, Department of Biophysics\n"
"; ----------------------------------------------\n"
"; ATTENTION! Do not forget to use the proper version\n"
"; of the force field fork (not less than revision). \n"
"; Watch for corresponding force field at: \n"
";            bitbucket.com/comcon1\n"
"; ----------------------------------------------\n"
"; Please ascertain that the topology is valid. We \n"
"; do not guarantee that. If you find that something\n"
"; is wrong, please report us to " PACKAGE_BUGREPORT "\n";

extern void save_topology(t_topology &, const char *) throw (t_exception);

extern void save_topology_rtp(t_topology &, const char *) throw (t_exception);

extern void load_topology(t_topology &, const char *) throw (t_exception);

extern void load_lack(t_topology &, const char *) throw (t_exception);

extern void check_topology(t_topology &) throw (t_exception);

extern void load_struct_stream(t_topology &, t_iformat, std::istream *) throw (t_exception);
extern void load_struct_fname(t_topology &, t_iformat, const char *) throw (t_exception);

extern void save_struct(t_topology &, t_oformat, const char *) throw (t_exception);

extern void save_lack(t_topology &, const char *);

#if ENABLE_GAMESS_FEATURES
extern void load_hessian(ublas::matrix<double>&, const char *) throw (t_exception);
#endif

//TODO: also serialization should be included

}
#endif

