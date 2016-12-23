AUTOMAKE_OPTIONS = subdir-objects
AM_CPPFLAGS = -I../include
AM_CXXFLAGS = -DBOOST_LOG_DYN_LINK -std=c++11
bin_PROGRAMS = programs/tpprenum programs/tppmktop
if GAMESSFEATURES
bin_PROGRAMS += programs/tppreadbonds
programs_tppreadbonds_SOURCES = programs/tppreadbonds.cpp structio.cpp topio.cpp hessian_accept.cpp calc.cpp runtime.cpp  logger.cpp strutil.cpp   db_base.cpp atom_definer.cpp bond_definer.cpp bond_matrix.cpp core.cpp
programs_tppreadbonds_LDADD = -lpthread -lboost_system -lboost_thread  -lboost_log_setup  -lboost_log  -llapack -lmysqlpp
endif
if SELFTESTFEATURES
bin_PROGRAMS += programs/tppmktest
programs_tppmktest_SOURCES = programs/tppmktest.cpp testcase.cpp runtime.cpp structio.cpp topio.cpp strutil.cpp tppnames.cpp logger.cpp db_base.cpp atom_definer.cpp bond_definer.cpp bond_matrix.cpp core.cpp
programs_tppmktest_LDADD = -lpthread -lboost_system -lboost_thread  -lboost_log_setup  -lboost_log  -lmysqlpp
endif
programs_tpprenum_SOURCES = programs/tpprenum.cpp strutil.cpp structio.cpp tppnames.cpp pdbutils.cpp exceptions.cpp logger.cpp bond_matrix.cpp core.cpp
programs_tppmktop_SOURCES = programs/tppmktop.cpp runtime.cpp structio.cpp topio.cpp strutil.cpp tppnames.cpp exceptions.cpp db_base.cpp atom_definer.cpp bond_definer.cpp bond_matrix.cpp core.cpp
programs_tppmktop_LDADD = -lpthread -lboost_system -lboost_thread  -lboost_log_setup  -lboost_log -lmysqlpp
