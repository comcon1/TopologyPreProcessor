#include "mysql.h"
#include "boost/multi_index_container.hpp"
#include "db_scanner.hpp"
#include "openbabel/obconversion.h"
#include "openbabel/obiter.h"
#include "boost/lambda/lambda.hpp"
#include "boost/lambda/if.hpp"

//#define CDB

namespace tpp {
  using namespace OpenBabel;

// common place procedures...
bond_definer::bond_definer(t_input_params _par, t_topology &_tp) throw (t_exception): 
 tp(_tp), par(_par), con(new mysqlpp::Connection(false)){
   connect_db();
}

bool bond_definer::connect_db() throw (t_exception) {
     string ffname = PARAM_READ(par,"ffname");
     int cat, cbon, cang, cdih, cnb;
     con->connect(
         PARAM_READ(par,"dbname").c_str(),
         (PARAM_READ(par,"host")+string(":")+PARAM_READ(par,"port")).c_str(),
         PARAM_READ(par,"user").c_str(),
         PARAM_READ(par,"password").c_str()
         );
     // connection established
     if (!con->connected()) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::bond_definer::connect_db");
        PARAM_ADD(params, "error", "SQL connection error");
        PARAM_ADD(params, "sql_error", "Cann't connect to DB!" );
        throw t_sql_exception("SQL connection failed!", params);
     }
     // get ffid
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      qu << format("SELECT id FROM forcefield WHERE name='%1$s'") % ffname.c_str();
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
        throw t_exception("Force field not found!", params);
      }
      ffid = res.at(0)["id"];
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
void bond_definer::fill_bonds() throw (t_exception) {
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
SELECT f,c1,c2 \
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
        throw t_exception("Bond definition error!", params);
      }
      runtime.log_write(string("Bond not found between ") + typ1 + " and " + typ2 + "!\n");
    }
     t_top_coord   tpc;
     t_top_element tel;
     tel.defname = string("dfTPP_bon_") + lexical_cast<string>(co); 
     tel.i = it->GetBeginAtomIdx();
     tel.j = it->GetEndAtomIdx();
     tpc.defname = tel.defname;
     tpc.type = TPP_TTYPE_BON;
     tpc.f =  (res.num_rows() > 0) ? res.at(0)["f"]  : -1;
     tpc.c0 = (res.num_rows() > 0) ? res.at(0)["c1"] : 0.00;
     tpc.c1 = (res.num_rows() > 0) ? res.at(0)["c2"] : 0.00;
     tp.elements.push_back(tel);
     tp.parameters.insert(tpc);
  } // end FOR_BONDS_OF_MOL
  
  
  // clearing the same bond types  
  for (t_top_map::nth_index<1>::type::iterator it = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it) 
   if ( tp.elements.get<1>().find(it->defname) != tp.elements.get<1>().end() ) {
    for (t_top_map::nth_index<1>::type::iterator it0 = tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
        it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it0) 
      // if defines are equal ;-)
      if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && (it0->f == it->f) &&
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end() ) 
        ) { 
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

void bond_definer::fill_angles() throw (t_exception) {
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
SELECT f,c1,c2 \
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
        throw t_exception("Bond definition error!", params);
      }
      runtime.log_write(string("Angle between ") + typ1 + ", " + typ2 + ", " + typ3 + " not found!\n");
    }
     t_top_coord   tpc;
     t_top_element tel;
     tel.defname = string("dfTPP_ang_") + lexical_cast<string>(co); 
     tel.i = (*it)[1]+1;
     tel.j = (*it)[0]+1;
     tel.k = (*it)[2]+1;
     tpc.defname = tel.defname;
     tpc.type = TPP_TTYPE_ANG;
     tpc.f =  (res.num_rows() > 0) ? res.at(0)["f"] : -1;
     tpc.c0 = (res.num_rows() > 0) ? res.at(0)["c1"] : 0.00;
     tpc.c1 = (res.num_rows() > 0) ? res.at(0)["c2"] : 0.00;
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
      if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && (it0->f == it->f) &&
           (tp.elements.get<1>().find(it0->defname) !=  tp.elements.get<1>().end() ) && ( it->f != -1) 
         ) { 
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

void bond_definer::fill_dihedrals() throw (t_exception) {
  mysqlpp::Query qu = con->query();
  MYSQLPP_RESULT res;
  mysqlpp::Row    row;
  string query;
  int co = 0;
  // including single pair definition
  { 
    t_top_coord tpc;
    tpc.type = TPP_TTYPE_PAIR;
    tpc.defname = "ONE_PAIR";
    tpc.f = 1;
    tp.parameters.insert(tpc);
  }
  // end including single pair
  FOR_TORSIONS_OF_MOL(it,tp.mol) {
    co++;
    qu.reset();
    string typ1 = namemap.find( tp.atoms.find((*it)[0]+1)->atom_type) ->second; 
    string typ2 = namemap.find( tp.atoms.find((*it)[1]+1)->atom_type)->second;
    string typ3 = namemap.find( tp.atoms.find((*it)[2]+1)->atom_type) ->second; 
    string typ4 = namemap.find( tp.atoms.find((*it)[3]+1)->atom_type) ->second; 
    qu << format("\
SELECT f,c1,c2,c3,c4,c5,c6 \
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
    if (res.num_rows() == 0 ) {
//      if (PARAM_EXISTS(par, "noqalculate")) {
//        t_input_params params;
//        PARAM_ADD(params, "procname", "tpp::bond_definer::fill_angles");
//        PARAM_ADD(params, "error", string("Angle not found between ") + typ1 + " and " + typ2 + " and " + typ3 + "!" );
//        throw t_exception("Bond definition error!", params);
//      } 
//      Making pair instead of exception
      t_top_element tel;
      tel.defname = "ONE_PAIR";
      tel.i = (*it)[0]+1;
      tel.j = (*it)[3]+1;
      tp.elements.push_back(tel);
    }
     t_top_coord   tpc;
     t_top_element tel;
     tel.defname = string("dfTPP_dih_") + lexical_cast<string>(co); 
     tel.i = (*it)[0]+1;
     tel.j = (*it)[1]+1;
     tel.k = (*it)[2]+1;
     tel.l = (*it)[3]+1;
     tpc.defname = tel.defname;
     if (res.num_rows() > 0) {
     tpc.f = (int)res.at(0)["f"];
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
  }
     else {
       tpc.f = -1;
       tpc.c0 = 0.00;
       tpc.c1 = 0.00;
       tpc.c2 = 0.00;
       tpc.c3 = 0.00;
       tpc.c4 = 0.00;
     }
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
      if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && 
           (it0->c2 == it->c2) && (it0->c3 == it->c3) && (it0->c4 == it->c4) &&
           (it0->c5 == it->c5) && (it0->f == it->f) && 
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end())
          ) { 
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

void bond_definer::fill_special() throw (t_exception) {
}
void bond_definer::bond_align() throw (t_exception) {
  try {
   fill_bonds();
   fill_angles();
   fill_dihedrals();
  } catch (t_sql_exception e) { e.fix_log(); }
    catch (t_exception e) { e.fix_log(); }
}
void bond_definer::log_needed_bonds() {
}

}

