
#include "mysql.h"
#include "testcase.hpp"

namespace tpp {

      /*
       * Function writes a whole structure file as a testcase
       * for further run tests.
       *
       * -- read filename, integer charge and comment string
       */
      void db_testcase::write_testmolecule(char *_fn, int _ch, char *_com) throw (t_exception) {

        return;

      }

      /*
      * Function writes a whole topology as a testcase
      * for further checking of DB and TPPMKTOP corectness.
      *
      * -- read t_topology object to write into DB table `selftest_cases`
      */
      void db_testcase::write_testcase(const t_topology &_tp) throw (t_exception) {

        //TODO: implement

        return;
      }

      /*
       * Function writes the testmolecule into the file.
       *
       * -- read molid (record ID in DB) and filename to write into
       * -- returns pair<charge, comment>
       */
      std::pair<int, string> db_testcase::save_testmolecule(int _ch, char *_fn) throw (t_exception) {

        //TODO: implement

        std::pair<int, string> ans(0, "");
        return ans;
      }

      /*
       * Perform test and print detailed information on a test.
       *
       * -- read molid and t_topology object
       * -- return test result
       */
      bool db_testcase::test_topology(int _molid, const t_topology &_tp) {

        //TODO: implement
        
        return false;
      }


}
