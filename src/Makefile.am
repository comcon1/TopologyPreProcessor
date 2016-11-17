AUTOMAKE_OPTIONS = subdir-objects
AM_CPPFLAGS = -I../include
AM_CXXFLAGS = -std=c++11
bin_PROGRAMS = programs/tpprenum programs/tppmktop
if GAMESSFEATURES
bin_PROGRAMS += programs/tppreadbonds
programs_tppreadbonds_SOURCES = programs/tppreadbonds.cpp topio.cpp hessian_accept.cpp calc.cpp runtime.cpp db_scanner.cpp dbinfo.cpp strutil.cpp db_bonds.cpp smartfit.cpp 
programs_tppreadbonds_LIBADD = -llapack -lmysqlpp
endif
if SELFTESTFEATURES
bin_PROGRAMS += programs/tppmktest
programs_tppmktest_SOURCES = programs/tppmktest.cpp testcase.cpp runtime.cpp db_scanner.cpp dbinfo.cpp smartfit.cpp db_bonds.cpp topio.cpp strutil.cpp tppnames.cpp 
programs_tppmktest_LIBADD = -lmysqlpp
endif
programs_tpprenum_SOURCES = programs/tpprenum.cpp strutil.cpp topio.cpp runtime.cpp pdbutils.cpp exceptions.cpp
programs_tppmktop_SOURCES = programs/tppmktop.cpp runtime.cpp db_scanner.cpp dbinfo.cpp smartfit.cpp db_bonds.cpp topio.cpp strutil.cpp tppnames.cpp exceptions.cpp
programs_tppmktop_LIBADD = -lmysqlpp
