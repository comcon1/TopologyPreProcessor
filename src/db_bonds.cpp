#include "mysql.h"

#include "runtime.hpp"

#include "boost/multi_index_container.hpp"
#include "db_scanner.hpp"
#include "tppnames.hpp"
#include "openbabel/obconversion.h"
#include "openbabel/obiter.h"

//#define CDB

namespace tpp {
  using namespace OpenBabel;

  // common place procedures...
  bond_definer::bond_definer(t_input_params _par, t_topology &_tp) throw (Exception):
    db_base(_par), tp(_tp) {
      connect_db();
  }

  bond_definer::~bond_definer() {        
    qalcfile.close();
  }


bool bond_definer::connect_db() throw (Exception) {
     db_base::connect_db();

     string ffname = PARAM_READ(par,"ffname");
     int cat, cbon, cang, cdih, cnb;

     // get ffid
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      qu << format("SELECT id, generate_pairs FROM forcefield WHERE name='%1$s'") % ffname.c_str();
      res = qu.store();
      if (!res) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::bond_definer::connect_db");
        PARAM_ADD(params, "error", "SQL query error");
        PARAM_ADD(params, "sql_error", qu.error() );
        throw t_sql_exception("SQL query failed!", params);
      }
      if (res.num_rows() == 0) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::bond_definer::connect_db");
        PARAM_ADD(params, "error", "Error in parameters");
        throw Exception("Force field not found!", params);
      }
      ffid = res.at(0)["id"];
      genpairs = (bool) (res.at(0)["generate_pairs"]);
      if (genpairs && (PARAM_READ(cmdline, "verbose_flag") == "on") ) {
        cout << "1-4 pair generation is required for FF." << endl;
      }
      qu.reset();
      
      // molecule stuff ))
      string query;
      FOR_ATOMS_OF_MOL(it,tp.mol) {
        t_atom_array::iterator pa = tp.atoms.find(it->GetIdx());
        namemap.insert(pair<string, string>(
              pa->atom_type, string("") )
            );
      }
      query = string("('quququq'");
      for (map<string,string>::iterator ii = namemap.begin(); ii != namemap.end(); ++ii) {
        query += ( string(",'") + ii->first + "'" );
      }
      query = string("SELECT `uname`,`name` FROM `atoms` WHERE `ffield` = ") + lexical_cast<string>(ffid) 
          + " and `uname` IN " + query + ")";
      // mysql stuff
      qu << query;
      res = qu.store();
      if (!res) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::bond_definer::connect_db");
        PARAM_ADD(params, "error", "SQL query error");
        PARAM_ADD(params, "sql_error", qu.error() );
        throw t_sql_exception("SQL query failed!", params);
      }

      for(mysqlpp::Row::size_type co = 0; co < res.num_rows(); ++co) {
        row = res.at(co);
        BOOST_CHECK(namemap.find(string(row["uname"].c_str())) != namemap.end() );
        namemap.find(string(row["uname"].c_str()))->second = row["name"].c_str();
      }
      return true;
}

// special functions))
void bond_definer::fill_bonds() throw (Exception) {
  mysqlpp::Query qu = con->query();
  MYSQLPP_RESULT res;
  mysqlpp::Row    row;
  string query;
  int co = 0;
  runtime.log_write("Defining bonds:\n");
  FOR_BONDS_OF_MOL(it,tp.mol) {
    co++;
    qu.reset();
    runtime.log_write(lexical_cast<string>(it->GetBeginAtomIdx())+ "->"+lexical_cast<string>(it->GetEndAtomIdx())+"\n");
    string typ1 = namemap.find( tp.atoms.find(it->GetBeginAtomIdx())->atom_type )->second;
    string typ2 = namemap.find( tp.atoms.find(it->GetEndAtomIdx())->atom_type) ->second; 
    qu << format("\
SELECT id,f,c1,c2 \
FROM bonds \
WHERE (bonds.ffield = %3$d) AND \
  ( (bonds.i = '%1$s' and bonds.j = '%2$s') OR \
    (bonds.i = '%2$s' and bonds.j = '%1$s')\
  )") % typ1 % typ2 % ffid;
     res = qu.store();
     if (!res) {
       t_input_params params;
       PARAM_ADD(params, "procname", "tpp::bond_definer::fill_bonds");
       PARAM_ADD(params, "error", "SQL query error");
       PARAM_ADD(params, "sql_error", qu.error() );
       throw t_sql_exception("SQL query failed!", params);
     }
     if ( res.num_rows() == 0 ) {
       if (PARAM_EXISTS(par, "noqalculate")) {
         t_input_params params;
         PARAM_ADD(params, "procname", "tpp::bond_definer::fill_bonds");
         PARAM_ADD(params, "error", string("Bond not found between ") + typ1 + " and " + typ2 + "!" );
         throw Exception("Bond definition error!", params);
       } else {
         ostringstream os;
         os << format("Bond not found between %1$s and %2$s!") % typ1 % typ2;
         runtime.log_write(os.str());
         if (PARAM_READ(cmdline, "verbose_flag") == "on") {
           cout << "[LACK] " << os.str() << endl;
         }
       }
       runtime.log_write(string("Bond not found between ") + typ1 + " and " + typ2 + "!\n");
     }
     t_top_coord   tpc;
     t_top_element tel;
     tel.i = it->GetBeginAtomIdx();
     tel.j = it->GetEndAtomIdx();
     tpc.type = TPP_TTYPE_BON;
     // FIXME: need to process exception when num_rows > 1
     tpc.f =  (res.num_rows() > 0) ? res.at(0)["f"]  : -1;
     tpc.c0 = (res.num_rows() > 0) ? res.at(0)["c1"] : 0.00;
     tpc.c1 = (res.num_rows() > 0) ? res.at(0)["c2"] : 0.00;
     tpc.dbid = (res.num_rows() > 0) ? res.at(0)["id"] : -1;
     tel.defname = ttc_name_generator(tpc).set_btypes({typ1,typ2}).getName(); 
     tpc.defname = tel.defname;
     tp.elements.push_back(tel);
     tp.parameters.insert(tpc);
  } // end FOR_BONDS_OF_MOL
  
  
  // clearing the same bond types  
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it) 
   if ( tp.elements.get<1>().find(it->defname) != tp.elements.get<1>().end() ) {
    for (t_top_map::nth_index<1>::type::iterator it0 = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
        it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it0) 
      // if defines are equal numerically
/*
      if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && (it0->f == it->f) &&
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end() ) 
        ) */
      // if defines are equal in DB
      if ( (it0 != it) && (it0->dbid == it->dbid) && 
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end()) ) { 
        string lastdef = it0->defname;
        t_top_array::nth_index<1>::type::iterator fnd = tp.elements.get<1>().find(it0->defname);
        t_top_element emod = *fnd;
        emod.defname = it->defname;
        tp.elements.get<1>().replace(fnd, emod); // correctly replacement with multi_index
        }
   }
  // erasing non-meaning parameters 
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it) 
   if ( tp.elements.get<1>().find(it->defname) == tp.elements.get<1>().end() ) 
     tp.parameters.get<1>().erase(it);

}

void bond_definer::fill_angles() throw (Exception) {
  mysqlpp::Query qu = con->query();
  MYSQLPP_RESULT res;
  mysqlpp::Row    row;
  string query;
  int co = 0;
  FOR_ANGLES_OF_MOL(it,tp.mol) {
    co++;
    qu.reset();
    string typ2 = namemap.find( tp.atoms.find((*it)[0]+1)->atom_type)->second;
    string typ1 = namemap.find( tp.atoms.find((*it)[1]+1)->atom_type) ->second; 
    string typ3 = namemap.find( tp.atoms.find((*it)[2]+1)->atom_type) ->second; 
    qu << format("\
SELECT id,f,c1,c2 \
FROM angles \
WHERE (angles.ffield = %4$d) AND \
  ( (angles.i = '%1$s' and angles.j = '%2$s' and angles.k = '%3$s') OR \
    (angles.i = '%3$s' and angles.j = '%2$s' and angles.k = '%1$s')\
  )") % typ1 % typ2  % typ3 % ffid;
     res = qu.store();

     if (!res) {
       t_input_params params;
       PARAM_ADD(params, "procname", "tpp::bond_definer::fill_angles");
       PARAM_ADD(params, "error", "SQL query error");
       PARAM_ADD(params, "sql_error", qu.error() );
       throw t_sql_exception("SQL query failed!", params);
     }

     if (res.num_rows() == 0 ) {
 
       if (PARAM_EXISTS(par, "noqalculate")) {
         t_input_params params;
         PARAM_ADD(params, "procname", "tpp::bond_definer::fill_angles");
         PARAM_ADD(params, "error", string("Angle not found between ") + typ1 + " and " + typ2 + " and " + typ3 + "!" );
         throw Exception("Bond definition error!", params);
       } else {
          ostringstream os;
          os << format("Angle not found between %1$s, %2$s and %3$s!") % typ1 % typ2 % typ3;
          runtime.log_write(os.str());
          if (PARAM_READ(cmdline, "verbose_flag") == "on") {
            cout << "[LACK] " << os.str() << endl;
          }
       }
 
     }

     t_top_coord   tpc;
     t_top_element tel;
     tel.i = (*it)[1]+1;
     tel.j = (*it)[0]+1;
     tel.k = (*it)[2]+1;
     tpc.type = TPP_TTYPE_ANG;
     tpc.f    = (res.num_rows() > 0) ? res.at(0)["f"]  : -1;
     tpc.c0   = (res.num_rows() > 0) ? res.at(0)["c1"] : 0.00;
     tpc.c1   = (res.num_rows() > 0) ? res.at(0)["c2"] : 0.00;
     tpc.dbid = (res.num_rows() > 0) ? res.at(0)["id"] : -1;
     tel.defname = ttc_name_generator(tpc).set_btypes({typ1,typ2,typ3}).getName(); 
     tpc.defname = tel.defname;
     tp.elements.push_back(tel);
     tp.parameters.insert(tpc);

  } // end FOR_BONDS_OF_MOL

  // clearing the same bond types
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it) 
   if (tp.elements.get<1>().find(it->defname) !=  tp.elements.get<1>().end() ) {
    for (t_top_map::nth_index<1>::type::iterator it0 = tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
        it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it0) 
      // if defines are equal ;-)
      //if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && (it0->f == it->f) &&
      //     (tp.elements.get<1>().find(it0->defname) !=  tp.elements.get<1>().end() ) && ( it->f != -1) 
      //   ) {
      //  now see only at dbid 
      if ( (it0 != it) && (it0->dbid == it->dbid) && 
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end()) ) { 
        t_top_array::nth_index<1>::type::iterator fnd = tp.elements.get<1>().find(it0->defname);
        t_top_element emod = *fnd;
        emod.defname = it->defname;
        tp.elements.get<1>().replace(fnd, emod); // correctly replacement with multi_index
      }
  }

  // erasing non-meaning parameters 
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it) 
   if ( tp.elements.get<1>().find(it->defname) == tp.elements.get<1>().end() ) 
     tp.parameters.get<1>().erase(it);
}

void bond_definer::fill_special() throw (Exception) {
 // TODO: add special dihedrals according to SMARTS
}

void bond_definer::fill_dihedrals() throw (Exception) {
  mysqlpp::Query qu = con->query();
  MYSQLPP_RESULT res;
  mysqlpp::Row    row;
  string query;
  int co = 0;
  FOR_TORSIONS_OF_MOL(it,tp.mol) {
    // TODO: check if special dihedral doesn't match
    co++;
    qu.reset();
    string typ1 = namemap.find( tp.atoms.find((*it)[0]+1)->atom_type) ->second; 
    string typ2 = namemap.find( tp.atoms.find((*it)[1]+1)->atom_type)->second;
    string typ3 = namemap.find( tp.atoms.find((*it)[2]+1)->atom_type) ->second; 
    string typ4 = namemap.find( tp.atoms.find((*it)[3]+1)->atom_type) ->second; 
    qu << format("\
SELECT id,f,c1,c2,c3,c4,c5,c6 \
FROM dihedrals \
WHERE (dihedrals.ffield = %5$d) AND \
  ( (dihedrals.i IN ('%1$s','X') and dihedrals.j IN ('%2$s','X') and dihedrals.k IN ('%3$s','X') and dihedrals.l IN ('%4$s','X') ) OR \
    (dihedrals.i IN ('%4$s','X') and dihedrals.j IN ('%3$s','X') and dihedrals.k IN ('%2$s','X') and dihedrals.l IN ('%1$s','X') )\
  )") % typ1 % typ2  % typ3 % typ4 % ffid;
    res = qu.store();
    if (!res) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::bond_definer::fill_dihedrals");
      PARAM_ADD(params, "error", "SQL query error");
      PARAM_ADD(params, "sql_error", qu.error() );
      throw t_sql_exception("SQL query failed!", params);
    }
     t_top_coord   tpc;
     t_top_element tel;
     tel.i = (*it)[0]+1;
     tel.j = (*it)[1]+1;
     tel.k = (*it)[2]+1;
     tel.l = (*it)[3]+1;
     if (res.num_rows() > 0) {
        tpc.f = (int)res.at(0)["f"];
        tpc.dbid = (int)res.at(0)["id"];
        switch (tpc.f) {
          case 3: 
            tpc.type = TPP_TTYPE_RBDIH;
            tpc.c0 = (double)res.at(0)["c1"];
            tpc.c1 = (double)res.at(0)["c2"];
            tpc.c2 = (double)res.at(0)["c3"];
            tpc.c3 = (double)res.at(0)["c4"];
            tpc.c4 = (double)res.at(0)["c5"];
            tpc.c5 = (double)res.at(0)["c6"];
            break;
          case 2:
            tpc.type = TPP_TTYPE_IMPDIH;
            tpc.c0 = (double)res.at(0)["c1"];
            tpc.c1 = (double)res.at(0)["c2"];
            tpc.c2 = (double)res.at(0)["c3"];
            tpc.c3 = 0;
            tpc.c4 = 0;
            tpc.c5 = 0;
            break;
          case 1:
            tpc.type = TPP_TTYPE_SYMDIH;
            tpc.c0 = (double)res.at(0)["c1"];
            tpc.c1 = (double)res.at(0)["c2"];
            tpc.c2 = (double)res.at(0)["c3"];
            tpc.c3 = 0;
            tpc.c4 = 0;
            tpc.c5 = 0;
            break;
          default: BOOST_ERROR("Wrong dihedral type");
        };
     } else {
       tpc.f = -1;
       //TODO: tpc type of undefined dihedral should come from FF defaults
       tpc.type = TPP_TTYPE_RBDIH;
       tpc.dbid = -1;
       tpc.c0 = 0.00;
       tpc.c1 = 0.00;
       tpc.c2 = 0.00;
       tpc.c3 = 0.00;
       tpc.c4 = 0.00;


       if (PARAM_EXISTS(par, "noqalculate")) {
         t_top_element tel_;
         tel_.defname = "ONE_PAIR";
         tel_.i = (*it)[0]+1;
         tel_.j = (*it)[3]+1;
         tp.elements.push_back(tel_);
         
         cout << format("[LACK] Additional pair was added instead of dihedral: %1$d-%1$d-%1$d-%1$d. \n") % 
           tel.i % tel.j % tel.k % tel.l;

/* // make pair instead of error
 *
         t_input_params params;
         PARAM_ADD(params, "procname", "tpp::bond_definer::fill_dihedrals");
         ostringstream os;
         os << format("Angle not found between %1$s, %2$s, %3$s and %4$s!") % typ1 % typ2 % typ3 % typ4;
         PARAM_ADD(params, "error", os.str() );
         throw Exception("Bond definition error!", params);
*/
       } else {
         ostringstream os;
         os << format("Dihedral not found between %1$s, %2$s, %3$s and %4$s!") % typ1 % typ2 % typ3 % typ4;
         runtime.log_write(os.str());
         if (PARAM_READ(cmdline, "verbose_flag") == "on") {
           cout << "[LACK] " << os.str() << endl;
         }
       }


     }
     tel.defname = ttc_name_generator(tpc).set_btypes({typ1,typ2,typ3,typ4}).getName(); 
     tpc.defname = tel.defname;
     tp.elements.push_back(tel);
     tp.parameters.insert(tpc);
  } // end FOR_BONDS_OF_MOL

  // clearing the same bond types
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it) 
   if (tp.elements.get<1>().find(it->defname) != tp.elements.get<1>().end()) {
    for (t_top_map::nth_index<1>::type::iterator it0 = tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
        it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it0) 
      // if defines are equal ;-)
      // if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && 
      //     (it0->c2 == it->c2) && (it0->c3 == it->c3) && (it0->c4 == it->c4) &&
      //     (it0->c5 == it->c5) && (it0->f == it->f) && 
      //     (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end())
      //    ) { 
      if ( (it0 != it) && (it0->dbid == it->dbid) && 
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end()) ) { 
        t_top_array::nth_index<1>::type::iterator fnd = tp.elements.get<1>().find(it0->defname);
        t_top_element emod = *fnd;
        emod.defname = it->defname;
        tp.elements.get<1>().replace(fnd, emod); // correctly replacement with multi_index
      }
  }

  // erasing non-meaning parameters 
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it) 
   if ( tp.elements.get<1>().find(it->defname) == tp.elements.get<1>().end() ) 
     tp.parameters.get<1>().erase(it);

}

// use for clever choosing and posing improper dihedrals
void bond_definer::fill_impropers() throw (Exception,DbException) {
    runtime.log_write("Starting curious SMART-improper-dihedral fitting.\n");
    mysqlpp::Query qu = con->query();
    qu <<  format("SELECT ip.id, ip.PAT, ip.order, ia.name, ip.impid, \
          ia.f, ia.c1, ia.c2, ia.c3, ia.c4, ia.c5, ia.c6 \
        FROM improper_patterns as ip \
        RIGHT JOIN impropers as ia ON ia.id = ip.impid \
        WHERE ip.ffield = %1$d and ip.override = 0 ") % this->ffid;
    runtime.log_write("Loading IMPROPER patterns from DB..");
    cout << "IMPROPER patterns are loading. Please wait.." << flush;          
    MYSQLPP_RESULT res;
    res = qu.store();
    if (!res) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::bond_definer::fill_impropers");
      PARAM_ADD(params, "error", "SQL query error");
      PARAM_ADD(params, "sql_error", qu.error() );
      PARAM_ADD(params, "sql_query", qu.str() ); 
      throw DbException("SQL query failed!", params);
    }
    runtime.log_write("OK!\n");
    cout << " finished.\n" <<      
      "Starting SMART-fit." << endl;
    cout << ( format("Patterns checked: %1$4d.") % 0 ) << flush;
    std::ostringstream os;
    mysqlpp::Row::size_type co;
    mysqlpp::Row    row;
    OBSmartsPattern pat;
    vector<vector<int> > maplist;
    for(co=0; co < res.num_rows(); ++co) {
      os.str(""); os.clear();
      row = res.at(co);
      pat.Init(row["PAT"]);
      os << format("[OB] Process PAT: %1$s having %2$d atoms.\n") 
        % row["PAT"] % pat.NumAtoms();
      runtime.log_write(os.str());
      pat.Match(tp.mol);
      BOOST_CHECK(pat.NumAtoms() == 4);
      maplist.clear();
      maplist = pat.GetUMapList();
#ifdef CDB
        cout << "============> " << row["PAT"] << "\n"
        << "Matches: " << maplist.size() << endl;
#endif
      os << format("[OB] Pattern %1$s matches %2$d times.\n")
          % row["PAT"] % maplist.size();
      runtime.log_write(os.str());
      // manipulating matches
      for(int i=0;i<maplist.size();++i) {
              int oo = (int) (row["order"]);
              BOOST_CHECK(oo <= 4321 and oo >= 1234);
              t_top_element tel;
              BOOST_CHECK(oo % 10 <= 4); // 4312 -> 2
              tel.l = maplist[i][ oo % 10 - 1 ];
              oo = oo / 10;
              BOOST_CHECK(oo % 10 <= 4); // 431 -> 1
              tel.k = maplist[i][ oo % 10 - 1 ];
              oo = oo / 10;
              BOOST_CHECK(oo % 10 <= 4); // 43 -> 3
              tel.j = maplist[i][ oo % 10 - 1 ];
              oo = oo / 10;
              BOOST_CHECK(oo <= 4); // 4 -> 4
              tel.i = maplist[i][ oo - 1 ];
              tel.defname = (string) row["name"];
              BOOST_CHECK( 
                  (tp.atoms.find(tel.i) != tp.atoms.end()) &&
                  (tp.atoms.find(tel.j) != tp.atoms.end()) &&
                  (tp.atoms.find(tel.k) != tp.atoms.end()) &&
                  (tp.atoms.find(tel.l) != tp.atoms.end())
                );
#ifdef CDB
              cout << format(" --- %1$d  %2$d  %3$d  %4$d ---") % tel.i % tel.j % tel.k % tel.l << endl;
#endif
              // element is adding unconditionally
              tp.elements.push_back(tel);
              // top-coord is added only once
              if ( !tp.parameters.count(tel.defname) ) {
                  t_top_coord tpc;
                  tpc.type = TPP_TTYPE_SPECIMP;
                  tpc.defname = tel.defname;
                  tpc.f = (int) row["f"];
                  tpc.c0 = (double) row["c1"];
                  tpc.c1 = (double) row["c2"];
                  tpc.c2 = (double) row["c3"];
                  tpc.c3 = (double) row["c4"];
                  tpc.c4 = (double) row["c5"];
                  tpc.c5 = (double) row["c6"];
                  tp.parameters.insert(tpc);
              } // end if
      } // end for (matches)
#ifdef CDB
          cout << endl;
#endif
      cout << ( format("\b\b\b\b\b%1$4d.") % (int)co ) << flush;
    } // end for (rows) 

    cout << endl;

    //TODO: remove impropers if dublicate for specials

} // end fill_special

/*
 * Function makes special pairs and common 1-4 pairs - ALSO.
 */
void bond_definer::fill_pairs() throw (Exception) {
  if (genpairs) {
    // including single pair definition
        t_top_coord tpc;
        tpc.type = TPP_TTYPE_PAIR;
        tpc.defname = "ONE_PAIR";
        tpc.f = 1;
        tp.parameters.insert(tpc);
    // generating pairs
    cout << "Generating 1-4 pairs for FF needs.." << flush;
    FOR_TORSIONS_OF_MOL(it,tp.mol) {
      t_top_element tel;
      tel.defname = "ONE_PAIR";
      BOOST_CHECK(tp.atoms.find((*it)[0]+1) != tp.atoms.end());
      tel.i = tp.atoms.find((*it)[0]+1)->index; 
      BOOST_CHECK(tp.atoms.find((*it)[1]+1) != tp.atoms.end());
      tel.j = tp.atoms.find((*it)[3]+1)->index;
      tp.elements.push_back(tel);
    }
    cout << "ok." << endl;
  }
}

  void bond_definer::bond_align() throw (Exception) {
    try {
      fill_bonds();
      fill_angles();
      fill_special();
      fill_dihedrals();
      fill_impropers();
      fill_pairs();
    } 
    catch (DbException &e) { e.fix_log(); cout << "..something fails." << endl; }
    catch (t_sql_exception &e) { e.fix_log(); cout << "..something fails." << endl; }
    catch (Exception &e) { e.fix_log();  cout << "..something fails." << endl; }
  }

  void bond_definer::log_needed_bonds() {
  }

}

