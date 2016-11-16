//#include "mysql.h"
#include "testcase.hpp"
#include "sqltypes.h"

namespace tpp {


      /*
       * Function writes a whole structure file as a testcase
       * for further run tests.
       *
       * -- read filename, integer charge and comment string
       */
      int db_testcase::write_testmolecule(const char *_fn, int _ch, const char *_com) throw (t_exception) {

        // read file to string

        std::ifstream fs(_fn);
        BOOST_REQUIRE(fs.is_open());
        std::string pdbtext;
        char *cur = new char[1000];
        do {
          fs.getline(cur, 1000);
          pdbtext += cur;
          pdbtext += "\n";
        } while (fs.peek() != EOF); 
        delete[] cur;

        // write the record

        mysqlpp::Connection con("tppforcefield", "localhost", "tppuser", "estatic");
        mysqlpp::Query query = con.query();

        selftest row(0, this->ffid, pdbtext, _ch, _com); 
        query.insert(row);
        mysqlpp::SimpleResult res = query.execute();
        
        return res.insert_id();

      }

      /*
      * Function writes a whole topology as a testcase
      * for further checking of DB and TPPMKTOP corectness.
      *
      * -- read t_topology object to write into DB table `selftest_cases`
      */
      void db_testcase::write_testcase(const t_topology &_tp, int molid) throw (t_exception) {

        mysqlpp::Connection con("tppforcefield", "localhost", "tppuser", "estatic");
        mysqlpp::Query query = con.query();

        // working with atoms
        
        {
          std::vector<tpp::selftestcase_Atom> sca_ar;

          for (t_atom_array::iterator it = _tp.atoms.begin(); it != _tp.atoms.end(); ++it) {
            sca_ar.push_back( selftestcase_Atom(molid, this->ffid, *it) );
          }

          query.insert(sca_ar.begin(), sca_ar.end());
          mysqlpp::SimpleResult res = query.execute();
          cout << format("[TC] Atom records inserted: %1$d.") % res.rows() << endl;
        }

        // working with bonds

        if (_tp.parameters.get<1>().find(TPP_TTYPE_BON) != _tp.parameters.get<1>().end()) {

          std::vector<tpp::selftestcase_Bond> sca_br;

          for (TopMap::nth_index_iterator<1>::type it = _tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
              it != _tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it)
            for (TopArray::nth_index_iterator<1>::type it0 = _tp.elements.get<1>().lower_bound(it->defname);
                it0 != _tp.elements.get<1>().upper_bound(it->defname); ++it0) {
              sca_br.push_back( selftestcase_Bond(molid, this->ffid, *it, *it0) );
            }
          query.insert(sca_br.begin(), sca_br.end());
          mysqlpp::SimpleResult res = query.execute();
          cout << format("[TC] Bond records inserted: %1$d.") % res.rows() << endl;
        } else {
          cout << "[TC] The topology does not have any bonds (WHY?)." << endl;
        }

        // working with valent angles


        if (_tp.parameters.get<1>().find(TPP_TTYPE_ANG) != _tp.parameters.get<1>().end()) {

          std::vector<tpp::selftestcase_Angle> sca_gr;

          for (TopMap::nth_index_iterator<1>::type it = _tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
              it != _tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it)
            for (TopArray::nth_index_iterator<1>::type it0 = _tp.elements.get<1>().lower_bound(it->defname);
                it0 != _tp.elements.get<1>().upper_bound(it->defname); ++it0) {
                sca_gr.push_back( selftestcase_Angle(molid, this->ffid, *it, *it0) );
            }
          query.insert(sca_gr.begin(), sca_gr.end());
          mysqlpp::SimpleResult res = query.execute();
          cout << format("[TC] Angle records inserted: %1$d.") % res.rows() << endl;

        } else {
          cout << "[TC] The topology does not have any valent angles (WHY?)." << endl;
        }

        // working with proper dihedrals
      
        if ( ( _tp.parameters.get<1>().find(TPP_TTYPE_RBDIH) != _tp.parameters.get<1>().end() ) ||
            ( _tp.parameters.get<1>().find(TPP_TTYPE_SYMDIH) != _tp.parameters.get<1>().end() ) ||
            ( _tp.parameters.get<1>().find(TPP_TTYPE_IMPDIH) != _tp.parameters.get<1>().end() ) 
            ) {

          std::vector<tpp::selftestcase_Dihedral> sca_dr;

          for (TopMap::nth_index_iterator<1>::type it = _tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
              it != _tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it)
            for (TopArray::nth_index_iterator<1>::type it0 = _tp.elements.get<1>().lower_bound(it->defname);
                it0 != _tp.elements.get<1>().upper_bound(it->defname); ++it0) {
                sca_dr.push_back( selftestcase_Dihedral(molid, this->ffid, *it, *it0) );
            }
        } else {
          cout << "[TC] The topology does not have any proper dihedrals (WHY?)." << endl;
        }

        // working with improper dihedrals

        if ( _tp.parameters.get<1>().find(TPP_TTYPE_SPECIMP) != _tp.parameters.get<1>().end() ) {

          cout << "[TC] WARNING: IMPROPERS HAVE NOT BEEN IMPLEMENTED YET!" << endl;

          // this is a stub

         /* 
          std::vector<tpp::selftestcase_Dihedral> sca_ir;

          for (t_top_map::nth_index_iterator<1>::type it = _tp.parameters.get<1>().lower_bound(TPP_TTYPE_SPECIMP);
              it != _tp.parameters.get<1>().upper_bound(TPP_TTYPE_SPECIMP); ++it)
            for (t_top_array::nth_index_iterator<1>::type it0 = _tp.elements.get<1>().lower_bound(it->defname);
                it0 != _tp.elements.get<1>().upper_bound(it->defname); ++it0) {
                sca_ir.push_back( selftestcase_Improper(molid, this->ffid, *it, *it0) );
            }
          */
        } else {
          cout << "[TC] The topology does not have any improper dihedrals (it can be normal)." << endl;
        }
  
        return;
      }

      /*
       * Function writes the testmolecule into the file.
       *
       * -- read molid (record ID in DB) and filename to write into
       * -- returns pair<charge, comment>
       */
      std::pair<int, std::string> db_testcase::save_testmolecule(int _ch, const char *_fn) throw (t_exception) {

        //TODO: implement
        
        std::pair<int, std::string> ans(0, "");
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
