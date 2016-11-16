#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "core.hpp"
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

void save_topology(Topology &, const char *) ;

void save_topology_rtp(Topology &, const char *);

void load_topology(Topology &, const char *);

void load_lack(Topology &, const char *);

void check_topology(Topology &) ;

void load_struct_stream(Topology &, InputFormat, std::istream *) ;
void load_struct_fname(Topology &, InputFormat, const char *);

void save_struct(Topology &, OutputFormat, const char *) ;

void save_lack(Topology &, const char *);

#if ENABLE_GAMESS_FEATURES
extern void load_hessian(ublas::matrix<double>&, const char *);
#endif

//TODO: also serialization should be included

}
#endif

