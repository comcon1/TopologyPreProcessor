#ifndef TPP_TESTCASE_H
#define TPP_TESTCASE_H

#include "global.hpp"
#include "db_scanner.hpp"
#include <mysql++/mysql++.h>
#include <mysql++/ssqls.h>

namespace tpp {

  sql_create_4(selftest,
      1,4,
      mysqlpp::sql_int,     id,
      mysqlpp::sql_text,    pdb,
      mysqlpp::sql_int,     charge,
      mysqlpp::sql_varchar, comment);

  sql_create_15(selftestcase,
      1,15,
      mysqlpp::sql_int,        id,
      mysqlpp::sql_int,        molid,
      mysqlpp::sql_int,        ffid,
      mysqlpp::sql_enum,       type1,
      mysqlpp::sql_varchar,    type2,
      mysqlpp::sql_int,        i,
      mysqlpp::sql_int_null,   j,
      mysqlpp::sql_int_null,   k,
      mysqlpp::sql_int_null,   l,
      mysqlpp::sql_float,      c1,
      mysqlpp::sql_float_null, c2,
      mysqlpp::sql_float_null, c3,
      mysqlpp::sql_float_null, c4,
      mysqlpp::sql_float_null, c5,
      mysqlpp::sql_float_null, c6);


  class selftestcase_Atom : selftestcase {
    public:
      selftestcase_Atom(const mysqlpp::sql_int &molid, const mysqlpp::sql_int &ffid, const t_atom& a) :
        selftestcase( 0, molid, ffid, 
            "atom", a.atom_type, a.index, mysqlpp::null, mysqlpp::null, mysqlpp::null,
            a.charge, a.mass, mysqlpp::null, mysqlpp::null, mysqlpp::null, mysqlpp::null ) {

      }
  };


  class db_testcase : public db_base {
    public:
     
      /*
       * Function writes a whole structure file as a testcase
       * for further run tests.
       *
       * -- read filename, integer charge and comment string
       */
      void write_testmolecule(char *, int, char*) throw (t_exception);

      /*
      * Function writes a whole topology as a testcase
      * for further checking of DB and TPPMKTOP corectness.
      *
      * -- read t_topology object to write into DB table `selftest_cases`
      */
      void write_testcase(const t_topology&) throw (t_exception);

      /*
       * Function writes the testmolecule into the file.
       *
       * -- read molid (record ID in DB) and filename to write into
       * -- returns pair<charge, comment>
       */
      std::pair<int, string> save_testmolecule(int, char*) throw (t_exception);

      /*
       * Perform test and print detailed information on a test.
       *
       * -- read molid and t_topology object
       * -- return test result
       */
      bool test_topology(int, const t_topology&);

  };

}

#endif
