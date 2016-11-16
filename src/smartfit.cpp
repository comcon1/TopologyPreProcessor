#include "db_scanner.hpp"

#include "exceptions.hpp"
#include "runtime.hpp"

#include <openbabel/obiter.h>
#include <openbabel/parsmart.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

//#define CDB

using std::map;
using std::cout;
using std::flush;
using std::endl;
using std::string;
using std::vector;

using boost::format;
using boost::lexical_cast;

namespace tpp {

  using namespace OpenBabel;

void atom_definer::smart_cgnr() throw (Exception) {
      // zero all charge groups
      for (AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
          Atom nat = *it;
          nat.c_gnr = 0;
          tp.atoms.replace(it,nat);
      }
      // algo body
      try {
          runtime.log_write("Starting curious SMART-charge-group fitting.\n");
          mysqlpp::Query qu = con->query();
          qu << format("SELECT id,PAT,flag FROM chargegroups \
                  WHERE ffield = %1$d and flag = 1") % this->ffid;
          runtime.log_write("Loading CGNR patterns from DB..");
          cout << "CHARGEGROUP patterns are loading. Please wait.." << flush;          
          MYSQLPP_RESULT res;
          res = qu.store();
          if (!res) {
              t_input_params params;
              PARAM_ADD(params, "procname", "tpp::atom_definer::smart_cgnr");
              PARAM_ADD(params, "error", "SQL query error");
              PARAM_ADD(params, "sql_error", qu.error() );
              throw t_sql_exception("SQL query failed!", params);
          }
          runtime.log_write("OK!\n");
          cout << " finished.\n" <<      
                  "Starting SMART-fit." << endl;

          OBSmartsPattern pat;
          mysqlpp::Row    row;
          mysqlpp::Row::size_type co;
          vector<vector<int> > maplist;
          set<int> atoms_suite;
          cout << ( format("Patterns checked: %1$4d.") % 0 ) << flush;
          std::map<int, vector<string> > sized_patterns;
          std::ostringstream os;
          for (co=0; co < res.num_rows(); ++co) {
              os.str(""); os.clear();
              row = res.at(co);
              pat.Init(row["PAT"]);
              int pna = pat.NumAtoms();
              os << format("[DB] Getting PATTERN: %1$s having %2$d atoms.\n") 
                  % row["PAT"] % pna;
              runtime.log_write(os.str());
              if (sized_patterns.count(pna) == 0) 
                  sized_patterns[pna] = vector<string>();
              sized_patterns[pna].push_back(string(row["PAT"]));
              cout << ( format("\b\b\b\b\b%1$4d.") % (int)co ) << flush;
          }
          // Applying patterns from small size to big size
          TPP_INDEX curCG = 1;
          for(std::map<int, vector<string> >::iterator it 
                  = sized_patterns.begin(); it != sized_patterns.end(); 
                  ++it) {
              os.str(""); os.clear();
              os << format("Aplying PATTERNS of size %1$d (%2$d total): \n") 
                  % it->first % it->second.size();
              runtime.log_write(os.str());
              for(vector<string>::iterator ci = it->second.begin();
                      ci != it->second.end(); ++ci) {
                  os.str(""); os.clear();
                  pat.Init(ci->c_str());
                  pat.Match(tp.mol);
                  maplist.clear();
                  maplist = pat.GetUMapList();

                  os << format("[OB] Pattern %1$s matches %2$d times.\n")
                      % (*ci) % maplist.size();
                  runtime.log_write(os.str());

                  for(int i=0;i<maplist.size();++i) {
                      for (int j=0; j<maplist[i].size(); ++j) {
                          AtomArray::iterator cur_it = tp.atoms.find((int) (maplist[i][j]));
#ifdef CDB
                          cout << maplist[i][j] << " " << flush;
#endif
                          BOOST_CHECK(cur_it != tp.atoms.end());
                          Atom cur0 = *cur_it;
                          cur0.c_gnr = curCG; 
                          tp.atoms.replace(cur_it, cur0);
                      }
                      curCG++;
#ifdef CDB
                      cout << endl;
#endif
                  } 

                  cout << '.' << flush;
              }
          }
          { // independent block of CGR renumbering
            cout << endl << "Renumbering CGNR according to human-readable style.." << flush;
            int current_cgr = 1;
            std::set<TPP_INDEX> done_atoms;
            std::set<TPP_INDEX> _tempset;
            for(AtomArray::iterator it = tp.atoms.begin();
                    it != tp.atoms.end(); ++it) {
                if (done_atoms.count(it->index)) continue;
                // for atom as a separate group
                if (it->c_gnr == 0) {
                    Atom cur0 = *it;
                    cur0.c_gnr = current_cgr;
                    done_atoms.insert(cur0.index);
                    tp.atoms.replace(it, cur0);
                }
                // for atom - part of a group
                else {
#ifdef CDB
                    cerr << it->c_gnr << " " << tp.atoms.get<2>().count(it->c_gnr) << endl;
#endif
                    TPP_INDEX oldcgnr = it->c_gnr;
                    _tempset.clear();
                    for (AtomArray::nth_index_iterator<2>::type cit = tp.atoms.get<2>().lower_bound(oldcgnr);
                        cit != tp.atoms.get<2>().upper_bound(oldcgnr); ++cit) {
                        if (done_atoms.count(cit->index)) continue;
                        _tempset.insert(cit->index);
                        // modifying c_gnr in a separate cycle !! (Index Policy
                        // Needs: see multi_index documentation.. )
                    }
                    for (set<TPP_INDEX>::iterator ii = _tempset.begin(); ii != _tempset.end(); ++ii) {
                        AtomArray::iterator cit = tp.atoms.find(*ii);
                        BOOST_CHECK( cit != tp.atoms.end() );
                        Atom cur0 = *cit;
                        cur0.c_gnr = current_cgr;
                        tp.atoms.replace(cit, cur0);
#ifdef CDB
                        cerr << ".";
#endif
                    }
                    done_atoms.insert(_tempset.begin(), _tempset.end());
                }
                current_cgr += 1;
            }
            cout << "finished." << endl;
          } // ending CGNR renumbering block

      } catch (Exception &et) {
          cout << "-- CATCH AT SMART_CGNR! --" << endl;
          throw et;
      }
}

void atom_definer::smart_fit() throw (Exception) {

      // make zero-scored copy of scores map
      map<int, map<int, int> > sf_scores (scores);
      for (auto &i : sf_scores)
        for (auto &j : i.second) 
          j.second = 0;
      
      // next work with copied sf_scores
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

      // process every pattern while reading DB rows
      for(co=0; co < res.num_rows(); ++co) {
        row = res.at(co);
        pat.Init(row["PAT"]);
        std::ostringstream os;
        os << format("[OB] Process PAT: %1$s having %2$d atoms.\n") 
            % row["PAT"] % pat.NumAtoms();
        runtime.log_write(os.str());
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

        for(set_it = atoms_suite.begin(); set_it != atoms_suite.end(); ++set_it) {
          score_it = sf_scores.find(*set_it);
          BOOST_CHECK( score_it != sf_scores.end() );
          score_subit = score_it->second.find( row["atom_ids"] );
          if ( score_subit == score_it->second.end() ) {
            t_input_params params;
            PARAM_ADD(params, "procname", "tpp::atom_definer::smartfit");
            PARAM_ADD(params, "error", string("SMART atom pattern #") + 
                lexical_cast<string>(row["apid"]) + " in DB is for invalid atom type!");
            throw Exception("SMARTS-DB ERROR!", params);
          }
          // fixing bug with summarizing equimatching smarts weights
          if ( (int)row["good"] > score_subit->second ) {
            score_subit->second = row["good"];
          }

        }
        cout << ( format("\b\b\b\b\b%1$4d.") % (int)co ) << flush;
      } // scanning over every smarts
      cout << "\n";

      // apply smart scores to main scores map
      for (auto &i : sf_scores) 
        for (auto &j : i.second) 
          scores[i.first][j.first] += TPP_SMART_COEFF * j.second;
} // end smartfit

}

