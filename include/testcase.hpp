#ifndef TPP_TESTCASE_H
#define TPP_TESTCASE_H

#include "db_base.hpp"

namespace tpp {

  class db_testcase : public db_base {
    public:

      int ffid;

      /*
       * Initializing DB-INFO class.
       */
      db_testcase(t_input_params p) throw (t_exception): db_base(p) {
        this->ffid = boost::lexical_cast<int>(PARAM_READ(p, "ffid"));
        this->connect_db();
      }

      /*
       * Function writes a whole structure file as a testcase
       * for further run tests.
       *
       * -- read filename, integer charge and comment string
       * -- return inserted ID
       */
      int write_testmolecule(const char *, int, const char*) throw (t_exception);

      /*
      * Function writes a whole topology as a testcase
      * for further checking of DB and TPPMKTOP corectness.
      *
      * -- read t_topology object to write into DB table `selftest_cases` and
      *        integer molid
      */
      void write_testcase(const t_topology&, int) throw (t_exception);

      /*
       * Function writes the testmolecule into the file.
       *
       * -- read molid (record ID in DB) and filename to write into
       * -- returns pair<charge, comment>
       */
      std::pair<int, std::string> save_testmolecule(int, const char*) throw (t_exception);

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

