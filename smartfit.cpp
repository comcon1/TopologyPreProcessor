#include "db_scanner.hpp"
#include "openbabel/obiter.h"
#include "openbabel/parsmart.h"

#define CDB

namespace tpp {

  using namespace OpenBabel;

void atom_definer::smart_fit() throw (t_exception) {
      runtime.log_write("Starting curious SMART-fitting procedure.\n");
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      mysqlpp::Row::size_type co;
      OBSmartsPattern pat;
      vector<vector<int> > maplist;
      set<int> atoms_suite;
      set<int>::iterator set_it;
      map<int, map<int,int> >::iterator score_it;
      map<int, int> score_submap;
      map<int, int>::iterator score_subit;
      qu << format("\
SELECT atom_patterns.id as apid, PAT, pos, atom_ids, atoms.znuc AS znuc, good \
FROM atom_patterns \
RIGHT JOIN atoms ON atoms.id = atom_patterns.atom_ids \
WHERE  (not atom_patterns.group = 1) and (atoms.ffield = %1$d)") % this->ffid;
      runtime.log_write("Loading patterns from database...");
      cout << "Patterns are loading. Please wait.." << flush;
      res = qu.store();
      if (!res) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::atom_definer::smartfit");
        PARAM_ADD(params, "error", "SQL query error");
        PARAM_ADD(params, "sql_error", qu.error() );
        throw t_sql_exception("SQL query failed!", params);
      }
      runtime.log_write("OK!\n");
      cout << " finished." << endl;
      cout << "Starting SMART-fit." << endl;
      cout << ( format("Patterns checked: %1$4d.") % 0 ) << flush;
      for(co=0; co < res.num_rows(); ++co) {
        row = res.at(co);
        pat.Init(row["PAT"]);
        pat.Match(tp.mol);
        maplist.clear();
        atoms_suite.clear();
        maplist = pat.GetMapList();
#ifdef CDB
        cout << "============> " << row["PAT"] << "\n"
        << "Matches: " << maplist.size() << endl;
#endif
        for(int i=0;i<maplist.size();++i) {
          BOOST_CHECK( maplist[i].size() >= (int)(row["pos"]) );
          atoms_suite.insert( maplist[i][row["pos"]-1] );
        } 
/*
        if ( (int)row["znuc"] == 1 ) {
          for(set_it = atoms_suite.begin(); set_it != atoms_suite.end(); ++set_it) {
            FOR_NBORS_OF_ATOM(it, (tp.mol.GetAtom(*set_it))) {
              if (it->GetAtomicNum() == 1) {
                score_it = scores.find(*set_it);
                BOOST_CHECK( score_it != scores.end() );
                score_subit = score_it->second.find(row["atom_ids"]);
                BOOST_CHECK( score_subit != score_it->second.end() );
                int old = score_subit->second;
                score_it->second.erase(score_subit);
                score_it->second.insert(pair<int,int>(row["atom_ids"],old+10*row["good"]));
              }
            } // over hydrogens cycle)
          }
        } else { */
          for(set_it = atoms_suite.begin(); set_it != atoms_suite.end(); ++set_it) {
            score_it = scores.find(*set_it);
//            cout << *set_it << "." << flush;
//            cout << row["atom_ids"] << "." << endl;
            BOOST_CHECK( score_it != scores.end() );
            score_subit = score_it->second.find( row["atom_ids"] );
            if ( score_subit == score_it->second.end() ) {
              t_input_params params;
              PARAM_ADD(params, "procname", "tpp::atom_definer::smartfit");
              PARAM_ADD(params, "error", string("SMART atom pattern #") + 
                  lexical_cast<string>(row["apid"]) + " in DB is for invalid atom type!");
              throw t_exception("SMARTS-DB ERROR!", params);
            }
            int old = score_subit->second;
            score_it->second.erase( score_subit );
            score_it->second.insert( pair<int,int>(row["atom_ids"],
                   old+TPP_SMART_COEFF*row["good"]) );
          }
/*        }*/
          cout << ( format("\b\b\b\b\b%1$4d.") % (int)co ) << flush;
      } // scanning over every smarts
      cout << "\n";
}

}

