#include "bond_definer.hpp"
#include "logger.hpp"
#include "tppnames.hpp"
#include "strutil.hpp"

#include <mysql.h>
#include <assert.h>
#include <boost/multi_index_container.hpp>

#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

//#define CDB

namespace tpp {

  using std::map;
  using std::cout;
  using std::pair;
  using std::flush;
  using std::endl;
  using std::string;
  using std::vector;
  using std::ostringstream;

  using boost::format;
  using boost::lexical_cast;

  using namespace OpenBabel;

  // common place procedures...
  BondDefiner::BondDefiner(const DbBase::Settings& s1,
                           const BondDefiner::Settings& s2,
                           Topology &_tp) : DbBase(s1),
                                           bondSettings(s2),
                                           tp(_tp) {
    cout << endl;
    TPPI << "== Starting BondDefiner ==";
    cout << endl;
    connectDB();
  }

  BondDefiner::~BondDefiner() {
    qalcfile.close();
  }

  /// special connection function
  bool BondDefiner::connectDB() {
    DbBase::connectDB();

    genPairs = strutil::split(tp.ffdefaults, " ")[2] == "yes";
    TPPD << format("Generating pairs flag: %d") % genPairs;
    if (genPairs ) {
      TPPD << "1-4 pair generation is required for FF.";
    }

    int cat, cbon, cang, cdih, cnb;
    mysqlpp::Query qu = con->query();
    ostringstream os;
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;

    // molecule stuff ))
    FOR_ATOMS_OF_MOL(it,tp.mol) {
      AtomArray::iterator pa = tp.atoms.find(it->GetIdx());
      namemap.insert(pair<string, string>(
              pa->atom_type, string("") )
      );
    }
    TPPD << "Request information about every atom used.";
    os.str("");
    os << format("SELECT `uname`,`name` FROM `atoms` WHERE `ffield` = %1$d and `uname` IN ") % bondSettings.ffID;
    os << "('Fvyxag8T*8dgw'"; // phrase that never meet in real table
    for (auto ii: namemap) {
      os << ",'" << ii.first << "'";
    }
    os << ")";
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    qu << os.str();
    res = qu.store();
    if (!res) {
      SqlException e("SQL query failed!");
      e.add("procname", "tpp::BondDefiner::connectDB");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error());
      throw e;
    }

    for (mysqlpp::Row::size_type co = 0; co < res.num_rows(); ++co) {
      row = res.at(co);
      assert(namemap.find(string(row["uname"].c_str())) != namemap.end());
      namemap.find(string(row["uname"].c_str()))->second =
          row["name"].c_str();
    }
    return true;
  } // connectDB

  // implement filling bonds
  void BondDefiner::fillBonds() {
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    ostringstream os;
    int co = 0;
    TPPD << "Defining bonds ...";
    FOR_BONDS_OF_MOL(it, tp.mol) {
      co++;
      qu.reset();
      TPPD << format("%1$d->%2$d") % it->GetBeginAtomIdx() % it->GetEndAtomIdx();
      string typ1 = namemap.find( tp.atoms.find(it->GetBeginAtomIdx())->atom_type )->second;
      string typ2 = namemap.find( tp.atoms.find(it->GetEndAtomIdx())->atom_type) ->second;
      os.str("");
      os << format(
        " SELECT id,f,c1,c2 "
        " FROM bonds "
        " WHERE (bonds.ffield = %3$d) AND "
        "  ( (bonds.i = '%1$s' and bonds.j = '%2$s') OR "
        "    (bonds.i = '%2$s' and bonds.j = '%1$s') "
        "  )") % typ1 % typ2 % bondSettings.ffID;
      #ifdef SQLDEBUG
      TPPD << os.str();
      #endif // SQLDEBUG
      qu << os.str();
      res = qu.store();
      if (!res) {
        SqlException e("SQL query failed!");
        e.add("procname", "tpp::BondDefiner::fill_bonds");
        e.add("error", "SQL query error");
        e.add("sql_error", qu.error() );
        throw e;
      }
      if ( res.num_rows() == 0 ) {
        if (bondSettings.finalize) { //TODO: change behaviour
          Exception e("Bond definition error!");
          e.add("procname", "tpp::BondDefiner::fill_bonds");
          e.add("error", string("Bond not found between ") + typ1 + " and " + typ2 + "!" );
          throw e;
        } else {
          ostringstream os_;
          os_ << format("Bond not found between %1$s and %2$s!") % typ1 % typ2;
          TPPD << os_.str();
          if (bondSettings.verbose) {
            cout << "[LACK] " << os_.str() << endl;
          }
        }
        TPPD << format("Bond not found between %1$s and %2$s!") % typ1 % typ2;
      }
      TopCoord tpc;
      TopElement tel;
      tel.i = it->GetBeginAtomIdx();
      tel.j = it->GetEndAtomIdx();
      tpc.type = TPP_TTYPE_BON;
      // FIXME: need to process exception when num_rows > 1
      tpc.f = (res.num_rows() > 0) ? res.at(0)["f"] : -1;
      tpc.c0 = (res.num_rows() > 0) ? res.at(0)["c1"] : 0.00;
      tpc.c1 = (res.num_rows() > 0) ? res.at(0)["c2"] : 0.00;
      tpc.dbid = (res.num_rows() > 0) ? res.at(0)["id"] : -1;
      tel.defname = TTCNameGenerator(tpc).set_btypes( {typ1,typ2}).getName();
      tpc.defname = tel.defname;
      tp.elements.push_back(tel);
      tp.parameters.insert(tpc);
    } // end FOR_BONDS_OF_MOL

    // clearing the same bond types
    for (TopMap::nth_index<1>::type::iterator it =
        tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
        it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it)
      if (tp.elements.get<1>().find(it->defname)
          != tp.elements.get<1>().end()) {
        for (TopMap::nth_index<1>::type::iterator it0 =
            tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
            it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON);
            ++it0)
          // if defines are equal numerically
          /*
           if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && (it0->f == it->f) &&
           (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end() )
           ) */
          // if defines are equal in DB
          if ((it0 != it) && (it0->dbid == it->dbid)
              && (tp.elements.get<1>().find(it0->defname)
                  != tp.elements.get<1>().end())) {
            string lastdef = it0->defname;
            TopArray::nth_index<1>::type::iterator fnd =
                tp.elements.get<1>().find(it0->defname);
            TopElement emod = *fnd;
            emod.defname = it->defname;
            tp.elements.get<1>().replace(fnd, emod); // correctly replacement with multi_index
          }
      }
    // end clearing the same bond types

    // erasing non-meaning parameters
    for (TopMap::nth_index<1>::type::iterator it =
        tp.parameters.get<1>().lower_bound(TPP_TTYPE_BON);
        it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_BON); ++it)
      if (tp.elements.get<1>().find(it->defname)
          == tp.elements.get<1>().end())
        tp.parameters.get<1>().erase(it);

  } // end fillBonds

void BondDefiner::fillAngles() {
  mysqlpp::Query qu = con->query();
  ostringstream os;
  QueryResult res;
  mysqlpp::Row row;
  string query;
  int co = 0;
  FOR_ANGLES_OF_MOL(it,tp.mol) {
    co++;
    qu.reset();
    string typ2 = namemap.find( tp.atoms.find((*it)[0]+1)->atom_type)->second;
    string typ1 = namemap.find( tp.atoms.find((*it)[1]+1)->atom_type)->second;
    string typ3 = namemap.find( tp.atoms.find((*it)[2]+1)->atom_type)->second;
    os.str("");
    os << format("\
  SELECT id,f,c1,c2 \
  FROM angles \
  WHERE (angles.ffield = %4$d) AND \
    ( (angles.i = '%1$s' and angles.j = '%2$s' and angles.k = '%3$s') OR \
      (angles.i = '%3$s' and angles.j = '%2$s' and angles.k = '%1$s')\
    )") % typ1 % typ2 % typ3 % bondSettings.ffID;
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif
    qu << os.str();
    res = qu.store();

    if (!res) {
      SqlException e("SQL query failed!");
      e.add("procname", "tpp::BondDefiner::fill_angles");
      e.add("error", "SQL query error");
      e.add("query", os.str() );
      e.add("sql_error", qu.error() );
      throw e;
    }

    if (res.num_rows() == 0 ) {

      if (bondSettings.finalize) { //TODO: change behaviour
        Exception e("Bond definition error!");
        e.add("procname", "tpp::BondDefiner::fill_angles");
        e.add("error", string("Angle not found between ") + typ1 + " and " + typ2 + " and " + typ3 + "!" );
        throw e;
      } else {
        ostringstream os;
        os << format("Angle not found between %1$s, %2$s and %3$s!") % typ1 % typ2 % typ3;
        TPPD << os.str();
        if (bondSettings.verbose) {
          cout << "[LACK] " << os.str() << endl;
        }
      }

    }

    TopCoord tpc;
    TopElement tel;
    tel.i = (*it)[1]+1;
    tel.j = (*it)[0]+1;
    tel.k = (*it)[2]+1;
    tpc.type = TPP_TTYPE_ANG;
    tpc.f = (res.num_rows() > 0) ? res.at(0)["f"] : -1;
    tpc.c0 = (res.num_rows() > 0) ? res.at(0)["c1"] : 0.00;
    tpc.c1 = (res.num_rows() > 0) ? res.at(0)["c2"] : 0.00;
    tpc.dbid = (res.num_rows() > 0) ? res.at(0)["id"] : -1;
    tel.defname = TTCNameGenerator(tpc).set_btypes( {typ1,typ2,typ3}).getName();
    tpc.defname = tel.defname;
    tp.elements.push_back(tel);
    tp.parameters.insert(tpc);

  } // end FOR_ANGLES_OF_MOL

  // clearing the same angle types
  for (TopMap::nth_index<1>::type::iterator it =
      tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it)
    if (tp.elements.get<1>().find(it->defname)
        != tp.elements.get<1>().end()) {
      for (TopMap::nth_index<1>::type::iterator it0 =
          tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
          it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG);
          ++it0)
        // if defines are equal ;-)
        //if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) && (it0->f == it->f) &&
        //     (tp.elements.get<1>().find(it0->defname) !=  tp.elements.get<1>().end() ) && ( it->f != -1)
        //   ) {
        //  now see only at dbid
        if ((it0 != it) && (it0->dbid == it->dbid)
            && (tp.elements.get<1>().find(it0->defname)
                != tp.elements.get<1>().end())) {
          TopArray::nth_index<1>::type::iterator fnd =
              tp.elements.get<1>().find(it0->defname);
          TopElement emod = *fnd;
          emod.defname = it->defname;
          tp.elements.get<1>().replace(fnd, emod); // correctly replacement with multi_index
        }
    }
  // end clearing the same angle types

  // erasing non-meaning parameters
  for (TopMap::nth_index<1>::type::iterator it =
      tp.parameters.get<1>().lower_bound(TPP_TTYPE_ANG);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_ANG); ++it)
    if (tp.elements.get<1>().find(it->defname)
        == tp.elements.get<1>().end())
      tp.parameters.get<1>().erase(it);
} // end fillAngles

void BondDefiner::fillSpecial() {
  // TODO: add special dihedrals according to SMARTS
}

void BondDefiner::fillDihedrals() {
  mysqlpp::Query qu = con->query();
  ostringstream os;
  QueryResult res;
  mysqlpp::Row row;
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
    os.str("");
    os << format("\
  SELECT id,f,c1,c2,c3,c4,c5,c6 \
  FROM dihedrals \
  WHERE (dihedrals.ffield = %5$d) AND \
    ( (dihedrals.i IN ('%1$s','X') and dihedrals.j IN ('%2$s','X') and dihedrals.k IN ('%3$s','X') and dihedrals.l IN ('%4$s','X') ) OR \
      (dihedrals.i IN ('%4$s','X') and dihedrals.j IN ('%3$s','X') and dihedrals.k IN ('%2$s','X') and dihedrals.l IN ('%1$s','X') )\
    )") % typ1 % typ2 % typ3 % typ4 % bondSettings.ffID;
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    qu << os.str();
    res = qu.store();
    if (!res) {
      SqlException e("SQL query failed!");
      e.add("procname", "tpp::BondDefiner::fill_dihedrals");
      e.add("error", "SQL query error");
      e.add("query", os.str() );
      e.add("sql_error", qu.error() );
      throw e;
    }
    TopCoord tpc;
    TopElement tel;
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
        default: throw std::runtime_error("Wrong dihedral type");
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

      if (bondSettings.finalize) { //TODO: change behaviour
        TopElement tel_;
        tel_.defname = "ONE_PAIR";
        tel_.i = (*it)[0]+1;
        tel_.j = (*it)[3]+1;
        tp.elements.push_back(tel_);

        cout << format("[LACK] Additional pair was added instead of dihedral: %1$d-%1$d-%1$d-%1$d. \n") %
        tel.i % tel.j % tel.k % tel.l;

        /* // make pair instead of error
         *
         Parameters params;
         params.add("procname", "tpp::BondDefiner::fill_dihedrals");
         ostringstream os;
         os << format("Angle not found between %1$s, %2$s, %3$s and %4$s!") % typ1 % typ2 % typ3 % typ4;
         params.add("error", os.str() );
         throw Exception("Bond definition error!", params);
         */
      } else {
        ostringstream os;
        os << format("Dihedral not found between %1$s, %2$s, %3$s and %4$s!") % typ1 % typ2 % typ3 % typ4;
        TPPD << os.str();
        if (bondSettings.verbose) {
          cout << "[LACK] " << os.str() << endl;
        }
      }

    }
    tel.defname = TTCNameGenerator(tpc).set_btypes( {typ1,typ2,typ3,typ4}).getName();
    tpc.defname = tel.defname;
    tp.elements.push_back(tel);
    tp.parameters.insert(tpc);
  } // end FOR_TORSIONS_OF_MOL

  // clearing the same bond types
  for (TopMap::nth_index<1>::type::iterator it =
      tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it)
    if (tp.elements.get<1>().find(it->defname)
        != tp.elements.get<1>().end()) {
      for (TopMap::nth_index<1>::type::iterator it0 =
          tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
          it0 != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH);
          ++it0)
        // if defines are equal ;-)
        // if ( (it0 != it) && (it0->c0 == it->c0) && (it0->c1 == it->c1) &&
        //     (it0->c2 == it->c2) && (it0->c3 == it->c3) && (it0->c4 == it->c4) &&
        //     (it0->c5 == it->c5) && (it0->f == it->f) &&
        //     (tp.elements.get<1>().find(it0->defname) != tp.elements.get<1>().end())
        //    ) {
        if ((it0 != it) && (it0->dbid == it->dbid)
            && (tp.elements.get<1>().find(it0->defname)
                != tp.elements.get<1>().end())) {
          TopArray::nth_index<1>::type::iterator fnd =
              tp.elements.get<1>().find(it0->defname);
          TopElement emod = *fnd;
          emod.defname = it->defname;
          tp.elements.get<1>().replace(fnd, emod); // correctly replacement with multi_index
        }
    }
  // end clearing the same bond types

  // erasing non-meaning parameters
  for (TopMap::nth_index<1>::type::iterator it =
      tp.parameters.get<1>().lower_bound(TPP_TTYPE_RBDIH);
      it != tp.parameters.get<1>().upper_bound(TPP_TTYPE_SYMDIH); ++it)
    if (tp.elements.get<1>().find(it->defname)
        == tp.elements.get<1>().end())
      tp.parameters.get<1>().erase(it);

} // end fillTorsions

//! use for clever choosing and posing improper dihedrals
void BondDefiner::fillImpropers() {
  TPPI << " ----> Performing automatic improper dihedral assignment..";
  mysqlpp::Query qu = con->query();
  ostringstream os;
  os << format(
          "SELECT ip.id, ip.PAT, ip.order, ia.name, ip.impid, \
          ia.f, ia.c1, ia.c2, ia.c3, ia.c4, ia.c5, ia.c6 \
        FROM improper_patterns as ip \
        RIGHT JOIN impropers as ia ON ia.id = ip.impid \
        WHERE ip.ffield = %1$d and ip.override = 0 ")
          % bondSettings.ffID;
  #ifdef SQLDEBUG
  TPPD << os.str();
  #endif // SQLDEBUG
  qu << os.str();
  TPPI << "  Loading IMPROPER patterns from DB.";
  QueryResult res;
  res = qu.store();
  if (!res) {
    DbException e("SQL query failed!");
    e.add("procname", "tpp::BondDefiner::fill_impropers");
    e.add("error", "SQL query error");
    e.add("sql_error", qu.error());
    e.add("sql_query", qu.str());
    throw e;
  }
  TPPI << "  Processing loaded SMART patterns.";
  if (!bondSettings.verbose)
    cout << (format("  Patterns checked: %4d.") % 0) << flush;
  mysqlpp::Row::size_type co;
  mysqlpp::Row row;
  OBSmartsPattern pat;
  vector<vector<int> > maplist;
  for (co = 0; co < res.num_rows(); ++co) {
    os.str("");
    os.clear();
    row = res.at(co);
    pat.Init(row["PAT"]);
    os
        << format("  ** [OB] Process PAT: %1$s having %2$d atoms.")
            % row["PAT"] % pat.NumAtoms();
    TPPD << os.str();
    pat.Match(tp.mol);
    assert(pat.NumAtoms() == 4);
    maplist.clear();
    maplist = pat.GetUMapList();
    #ifdef CDB
    cout << "============> " << row["PAT"] << "\n"
    << "Matches: " << maplist.size() << endl;
    #endif
    os
        << format("  ** [OB] Pattern %1$s matches %2$d times.")
            % row["PAT"] % maplist.size();
    TPPD << os.str();
    // manipulating matches
    for (int i = 0; i < maplist.size(); ++i) {
      int oo = (int) (row["order"]);
      assert(oo <= 4321 and oo >= 1234);
      TopElement tel;
      assert(oo % 10 <= 4); // 4312 -> 2
      tel.l = maplist[i][oo % 10 - 1];
      oo = oo / 10;
      assert(oo % 10 <= 4); // 431 -> 1
      tel.k = maplist[i][oo % 10 - 1];
      oo = oo / 10;
      assert(oo % 10 <= 4); // 43 -> 3
      tel.j = maplist[i][oo % 10 - 1];
      oo = oo / 10;
      assert(oo <= 4); // 4 -> 4
      tel.i = maplist[i][oo - 1];
      tel.defname = (string) row["name"];
      assert(
          (tp.atoms.find(tel.i) != tp.atoms.end())
              && (tp.atoms.find(tel.j) != tp.atoms.end())
              && (tp.atoms.find(tel.k) != tp.atoms.end())
              && (tp.atoms.find(tel.l) != tp.atoms.end()));
      #ifdef CDB
      cout << format(" --- %1$d  %2$d  %3$d  %4$d ---") % tel.i % tel.j % tel.k % tel.l << endl;
      #endif
      // element is adding unconditionally
      tp.elements.push_back(tel);
      // top-coord is added only once
      if (!tp.parameters.count(tel.defname)) {
        TopCoord tpc;
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
    if (!bondSettings.verbose)
      cout << (format("\b\b\b\b\b%1$4d.") % (int) co) << flush;
  } // end for (rows)

  if (!bondSettings.verbose)
    cout << endl;

  //TODO: remove impropers if duplicate for specials
  TPPI << "                                                          .. DONE. <----";

} // end fill_special

//! Function makes special pairs and common 1-4 pairs - ALSO.
void BondDefiner::fillPairs() {
  if (genPairs) {
    // including single pair definition
    TopCoord tpc;
    tpc.type = TPP_TTYPE_PAIR;
    tpc.defname = "ONE_PAIR";
    tpc.f = 1;
    tp.parameters.insert(tpc);
    // generating pairs
    TPPI << "Generating 1-4 pairs.";
    FOR_TORSIONS_OF_MOL(it,tp.mol) {
      TopElement tel;
      tel.defname = "ONE_PAIR";
      assert(tp.atoms.find((*it)[0]+1) != tp.atoms.end());
      tel.i = tp.atoms.find((*it)[0]+1)->index;
      assert(tp.atoms.find((*it)[1]+1) != tp.atoms.end());
      tel.j = tp.atoms.find((*it)[3]+1)->index;
      tp.elements.push_back(tel);
    }
  }
} // end fillPairs

void BondDefiner::bondAlign() {
  try {
    fillBonds();
    fillAngles();
    fillSpecial(); // THAT ORDER?
    fillDihedrals();
    fillImpropers();
    fillPairs();
  } catch (const DbException &e) {
    cout << "..something fails." << endl;
    TPPD << e.what();
  } catch (const SqlException &e) {;
    cout << "..something fails." << endl;
    TPPD << e.what();
  } catch (const Exception &e) {
    cout << "..something fails." << endl;
    TPPD << e.what();
  }
}

void BondDefiner::log_needed_bonds() {
}

} // of namepsace tpp

