#include "atom_definer.hpp"

#include <mysql.h>

#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <logger.hpp>

namespace tpp {

  using std::map;
  using std::set;
  using std::pair;
  using std::cout;
  using std::endl;
  using std::flush;
  using std::string;
  using std::ostringstream;
  using std::vector;

  using boost::format;
  using boost::lexical_cast;
  using namespace boost::multi_index;

  using namespace OpenBabel;
  using namespace tpp::detail;

  /**
    *  << WRITE THE ALGORITHM HERE >>
    *
    *
    *
    *
    */

  /// implement non-bond map filling
  void AtomDefiner::fillNB() {
    cout << "Comparing atom by elements.." << endl;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;
    set<string> tmpv;
    map<int, set<string> > znucset; // map atomic number -> set of valence atomtypes
    // znucset is set to { ATOMIC NUMBER -> EMPTY SET }
    // it contains only ATOMIC NUMBERs which occurs in the mol
    FOR_ATOMS_OF_MOL(it, tp.mol) {
      znucset.insert(pair<int, set<string> >(it->GetAtomicNum(), set<string>() ) );
    }
    // perform query for every ATOMIC NUMBER
    cout << "Atomic numbers: [" << flush;
    for (auto &cur: znucset) {
      qu = con->query();
      qu << format(
        "SELECT `atoms`.`name` as nm FROM `atoms` "
        "WHERE `znuc` = %1$d AND `ffield` = %2$d "
        "GROUP BY `name`") % cur.first % atomSettings.ffID;
      res = qu.store();
      if (!res) {
        Exception e("Empty atom list resulted.");
        e.add("query", qu.str());
        throw e;
      }
      for (co = 0; co < res.num_rows(); ++co) {
        row = res.at(co);
        cur.second.insert( std::string(row["nm"]) );
      }
      cout << "." << flush;
    }
    cout << format("]\n"
      "%1$d queries proceeded on database\n") % znucset.size();
    cout << flush;
    // fill map to every atom
    FOR_ATOMS_OF_MOL(it, tp.mol) {
      nbSuite.insert(pair<int,set<string> >(it->GetIdx(), znucset[it->GetAtomicNum()]));
    }
  } // end fillNB

  /// implement existing bond counting
  void AtomDefiner::fillBonds() {
    cout << "Comparing atom by bonds.." << endl;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;

    typedef Spec2_<string> NMBondSpec;   /// bond type specified by atom valence names
    typedef Spec2_<int>    AIBondSpec;   /// bond type specified by atom numbers only
    set<NMBondSpec> tmpv, tmpv1;
    int qry = 0;
    map<AIBondSpec, set<NMBondSpec> > znucset; /// map: AtomNumber-AtomNumber -> AtomName-AtomName

    // all bond types are included once
    FOR_BONDS_OF_MOL(it, tp.mol) {
      int i0 = it->GetBeginAtom()->GetIdx(), j0 = it->GetEndAtom()->GetIdx(); // i,j
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum();
      AIBondSpec S(a1 <= a2 ? a1 : a2, a1 <= a2 ? a2 : a1);
      znucset.insert( pair<AIBondSpec, set<NMBondSpec> > ( S, set<NMBondSpec>() ) );
    }

    // perform query at every Bond type
    cout << "Bonds: [" << flush;
    for (auto &cur: znucset) {
      int i = cur.first.first(), j = cur.first.second();
        qu = con->query();
        if (i == j)
          qu << format(
            " SELECT "
            "   `bonds`.`id` AS `bid`, `bonds`.`i` AS `i`, `bonds`.`j` AS `j` "
            " FROM `bonds` "
            " RIGHT JOIN `atoms` as `iatoms` ON `bonds`.`i` = `iatoms`.`name` "
            " RIGHT JOIN `atoms` as `jatoms` ON `bonds`.`j` = `jatoms`.`name` "
            " WHERE `bonds`.`ffield` = %1$d "
            "   AND `iatoms`.`znuc` = %2$d AND `jatoms`.`znuc` = %2$d "
            " GROUP BY `bid`"
            ) % atomSettings.ffID % i;
        else
          qu << format(
            " SELECT "
            "       `bonds`.`id` AS `bid`, `bonds`.`i` AS `i`, `bonds`.`j` AS `j` "
            "     FROM bonds "
            "     RIGHT JOIN `atoms` as `iatoms` ON `bonds`.`i` = `iatoms`.`name` AND `iatoms`.`ffield` = %1$d "
            "     RIGHT JOIN `atoms` as `jatoms` ON `bonds`.`j` = `jatoms`.`name` AND `iatoms`.`ffield` = %1$d "
            "     WHERE `bonds`.`ffield` = %1$d "
            "       AND `iatoms`.`znuc` = %2$d AND `jatoms`.`znuc` = %3$d "
            "     GROUP BY `bid` "
            "     "
            " UNION "
            "     "
            " SELECT "
            "       `bonds`.`id` AS `bid`, `bonds`.`j` AS `i`, `bonds`.`i` AS `j` "
            "     FROM bonds "
            "     RIGHT JOIN `atoms` as `iatoms` ON `bonds`.`i` = `iatoms`.`name` AND `iatoms`.`ffield` = %1$d "
            "     RIGHT JOIN `atoms` as `jatoms` ON `bonds`.`j` = `jatoms`.`name` AND `iatoms`.`ffield` = %1$d "
            "     WHERE `bonds`.`ffield` = %1$d "
            "       AND `iatoms`.`znuc` = %3$d AND `jatoms`.`znuc` = %2$d "
            "     GROUP BY `bid` "
            ) % atomSettings.ffID % i % j;
        res = qu.store();
        if (!res) {
          Exception e("Empty bond list resulted.");
          e.add("query", qu.str());
          throw e;
        }
        for (co = 0; co < res.num_rows(); ++co) {
          row = res.at(co);
          cur.second.insert( NMBondSpec(string(row["i"]),string(row["j"]) ) );
        }
        cout << "." << flush;
    } // end for znucset
    cout << format("]\n"
              "%1$d queries proceeded on database\n") % znucset.size();
    cout << flush;

    // associate znucset bondsets with every bond
    FOR_BONDS_OF_MOL(it, tp.mol) {
      Spec2<int> idxBnd(it->GetBeginAtom()->GetIdx(), it->GetEndAtom()->GetIdx()); // bond i->j
      int i0 = idxBnd.first(), j0 = idxBnd.second(); // i,j
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum();
      if ( znucset.count(AIBondSpec(a1,a2)) ) {
        bondSuite.insert( pair<Spec2<int>, set<NMBondSpec> > (
           idxBnd, znucset[AIBondSpec(a1,a2)] ) );
      } else if ( znucset.count(AIBondSpec(a2,a1) ) ) {
        // swap first and second in the set
        tmpv1.clear();
        for (auto ii: znucset[AIBondSpec(a1,a2)])
          tmpv1.insert(NMBondSpec(ii.second(),ii.first()));
        bondSuite.insert( pair<Spec2<int>, set<NMBondSpec> > (
           idxBnd, tmpv1 ) );
      } else {
        /// This situation can not occur. Just test that the code is correct.
        TPPE << "No AIBondSpec was prepared!";
        assert(0);
      }
    } // end FOR_BONDS_OF_MOL

  } // end fillBonds

  ///TODO: REFACTOR AS fillBonds
  void AtomDefiner::fillAngles() {
    return ; // not ready yet
/*
    cout << "Comparing atom by angles.." << endl;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;
    set<Spec3_<int> > tmpv, tmpv1;
    map<Spec3_<int>, set<Spec3_<int> > > znucset;
    short qry = 0;
    cout << "Angles: [" << flush;
    FOR_ANGLES_OF_MOL(it, tp.mol){
    Spec3<int> A( (*it)[1]+1, (*it)[0]+1, (*it)[2]+1); // one-stabled angle definition
    int i0 = A.first(), j0 = A.second(), k0 = A.third();
    int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
    a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
    a3 = tp.mol.GetAtom(k0)->GetAtomicNum();
  #ifdef CDB
    cout << (format("#%1$d-%2$d-%3$d,(%4$d-%5$d-%6$d)") % i0 %j0 % k0 % a1 % a2 % a3) << flush;
  #endif
    Spec3_<int> S(a1 <= a3 ? a1 : a3, a2, a1 <= a3 ? a3 : a1);
    if ( znucset.find(S) == znucset.end() ) {
      qry++;
      qu = con->query();
      if (a1 == a3)
      qu << format("\
  SELECT \n\
    angles.id AS bid, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.i AND atoms.ffield = %1$d) AS aid1, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.j AND atoms.ffield = %1$d) AS aid2, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.k AND atoms.ffield = %1$d) AS aid3, \n\
    0 AS inv\n\
  FROM angles\n\
  WHERE angles.ffield = %1$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.i AND atoms.ffield = %1$d) = %2$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.j AND atoms.ffield = %1$d) = %3$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.k AND atoms.ffield = %1$d) = %2$d\n\
  ") % atomSettings.ffID % a1 % a2;
      else
      qu << format("\
  (SELECT \n\
    angles.id AS bid, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.i AND atoms.ffield = %1$d) AS aid1, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.j AND atoms.ffield = %1$d) AS aid2, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.k AND atoms.ffield = %1$d) AS aid3, \n\
    0 AS inv\n\
  FROM angles\n\
  WHERE angles.ffield = %1$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.i AND atoms.ffield = %1$d) = %2$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.j AND atoms.ffield = %1$d) = %3$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.k AND atoms.ffield = %1$d) = %4$d)\n\
  UNION \n\
  (SELECT \n\
    angles.id AS bid, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.i AND atoms.ffield = %1$d) AS aid1, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.j AND atoms.ffield = %1$d) AS aid2, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = angles.k AND atoms.ffield = %1$d) AS aid3, \n\
    1 AS inv\n\
  FROM angles\n\
  WHERE angles.ffield = %1$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.i AND atoms.ffield = %1$d) = %4$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.j AND atoms.ffield = %1$d) = %3$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = angles.k AND atoms.ffield = %1$d) = %2$d)\n\
  ") % atomSettings.ffID % a1 % a2 % a3;
      res = qu.store();
      assert(res);
      tmpv.clear();
      for (co = 0; row = res.at(co); ++co) {
        tmpv.insert((int)row["inv"] == 0 ?
            Spec3_<int>(row["aid1"],row["aid2"],row["aid3"]) :
            Spec3_<int>(row["aid3"],row["aid2"],row["aid1"])
        );
      }
  #ifdef CDB
      cout << "\ninserted: " << S.first() << ":" << S.second() << ":" << S.third() << " elements: " <<
      std::distance(tmpv.begin(),tmpv.end()) << endl;
  #endif

      if (S.first() != a1) {
        tmpv1.clear();
        for (set<Spec3_<int> >::iterator iii = tmpv.begin(); iii != tmpv.end(); ++iii)
        tmpv1.insert(Spec3_<int>(iii->third(),iii->second(), iii->first()));
        znucset.insert(pair<Spec3_<int>, set<Spec3_<int> > >(S, tmpv1));
      } else {
        znucset.insert(pair<Spec3_<int>, set<Spec3_<int> > >(S, tmpv));
      }
    }  // end find in database
    else {
  #ifdef CDB
      cout << "FOUND: " << znucset.find(S)->first.first()
      << "-" << znucset.find(S)->first.second()
      << "-" << znucset.find(S)->first.third() << endl;
  #endif
      tmpv1.clear();
      tmpv1 = znucset.find(S)->second;
      if (a1 == S.first()) {
        tmpv = tmpv1;
      } else { // turning over set
        tmpv.clear();
        for (set<Spec3_<int> >::iterator iii = tmpv1.begin(); iii != tmpv1.end(); ++iii)
          tmpv.insert(Spec3_<int> (iii->third(), iii->second(), iii->first()));
      }
    } // end find in cash set
    angleSuite.insert(pair<Spec3<int>, set<Spec3_<int> > >(A,tmpv));
    cout << "." << flush;
  #ifdef CDB
    cout << std::distance(tmpv.begin(),tmpv.end());
  #endif
  } // end every angle
    cout << format("]\n\
  %1$d queries proceeded on database\n") % qry;
    cout << flush;
    */
  }

  ///TODO: REFACTOR AS fillBonds
  void AtomDefiner::fillDihs() {
    return; //TODO: not ready yet
/*
    cout << "Comparing atom by dihedrals.." << endl;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;
    set<Spec4_> tmpv, tmpv1;
    map<Spec4_, set<Spec4_> > znucset;
    short qry = 0;
    cout << "Dihedrals: [" << flush;
    FOR_TORSIONS_OF_MOL(it, tp.mol){
    Spec4 A( (*it)[1]+1, (*it)[0]+1, (*it)[2]+1, (*it)[3]+1); // one-stabled angle definition
    int i0 = A.first(), j0 = A.second(), k0 = A.third(), l0 = A.fourth();
    int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
    a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
    a3 = tp.mol.GetAtom(k0)->GetAtomicNum(),
    a4 = tp.mol.GetAtom(l0)->GetAtomicNum();
  #ifdef CDB
    cout << (format("#%1$d-%2$d-%3$d,(%4$d-%5$d-%6$d)") % i0 %j0 % k0 % a1 % a2 % a3) << flush;
  #endif
    Spec4_ S(a1 <= a4 ? a1 : a4, a1 <= a4 ? a2 : a3, a1 <= a4 ? a3 : a2, a1 <= a4 ? a4 : a1);
    if ( znucset.find(S) == znucset.end() ) {
      qry++;
      qu = con->query();
      if (a1 == a4)
      qu << format("\
  SELECT \n\
    dihedrals.id AS did, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.i AND atoms.ffield = %1$d) AS aid1, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.j AND atoms.ffield = %1$d) AS aid2, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.k AND atoms.ffield = %1$d) AS aid3, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.l AND atoms.ffield = %1$d) AS aid4, \n\
    0 AS inv\n\
  FROM dihedrals\n\
  WHERE dihedrals.ffield = %1$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.i AND atoms.ffield = %1$d) = %2$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.j AND atoms.ffield = %1$d) = %3$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.k AND atoms.ffield = %1$d) = %4$d\n\
  AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.l AND atoms.ffield = %1$d) = %2$d\n\
  ") % atomSettings.ffID % a1 % a2 % a3;
      else
      qu << format("\
  (SELECT \n\
    dihedrals.id AS bid, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.i AND atoms.ffield = %1$d) AS aid1, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.j AND atoms.ffield = %1$d) AS aid2, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.k AND atoms.ffield = %1$d) AS aid3, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.l AND atoms.ffield = %1$d) AS aid4, \n\
    0 AS inv\n\
  FROM dihedrals\n\
  WHERE dihedrals.ffield = %1$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.i AND atoms.ffield = %1$d) = %2$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.j AND atoms.ffield = %1$d) = %3$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.k AND atoms.ffield = %1$d) = %4$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.l AND atoms.ffield = %1$d) = %5$d)\n\
  UNION \n\
  (SELECT \n\
    dihedrals.id AS bid, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.i AND atoms.ffield = %1$d) AS aid1, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.j AND atoms.ffield = %1$d) AS aid2, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.k AND atoms.ffield = %1$d) AS aid3, \n\
    ( SELECT MIN( id ) FROM atoms WHERE atoms.name = dihedrals.l AND atoms.ffield = %1$d) AS aid4, \n\
    1 AS inv\n\
  FROM dihedrals\n\
  WHERE dihedrals.ffield = %1$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.i AND atoms.ffield = %1$d) = %5$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.j AND atoms.ffield = %1$d) = %4$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.k AND atoms.ffield = %1$d) = %3$d\n\
   AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = dihedrals.l AND atoms.ffield = %1$d) = %2$d)\n\
  ") % atomSettings.ffID % a1 % a2 % a3 % a4;
  //        cout << endl << qu.preview() << endl;
      res = qu.store();
      assert(res);
      tmpv.clear();
      for (co = 0; row = res.at(co); ++co) {
        tmpv.insert((int)row["inv"] == 0 ?
            Spec4_(row["aid1"],row["aid2"],row["aid3"],row["aid4"]) :
            Spec4_(row["aid4"],row["aid3"],row["aid2"],row["aid1"])
        );
      }
  #ifdef CDB
      cout << "\ninserted: " << S.first() << ":" << S.second() << ":" << S.third() << " elements: " <<
      std::distance(tmpv.begin(),tmpv.end()) << endl;
  #endif

      if (S.first() != a1) {
        tmpv1.clear();
        for (set<Spec4_>::iterator iii = tmpv.begin(); iii != tmpv.end(); ++iii)
        tmpv1.insert(Spec4_(iii->fourth(),iii->third(),iii->second(), iii->first()));
        znucset.insert(pair<Spec4_, set<Spec4_> >(S, tmpv1));
      } else {
        znucset.insert(pair<Spec4_, set<Spec4_> >(S, tmpv));
      }
    }  // end find in database
    else {
  #ifdef CDB
      cout << "FOUND: " << znucset.find(S)->first.first()
      << "-" << znucset.find(S)->first.second()
      << "-" << znucset.find(S)->first.third() << endl;
  #endif
      tmpv1.clear();
      tmpv1 = znucset.find(S)->second;
      if (a1 == S.first()) {
        tmpv = tmpv1;
      } else { // turning over set
        tmpv.clear();
        for (set<Spec4_>::iterator iii = tmpv1.begin(); iii != tmpv1.end(); ++iii)
        tmpv.insert(Spec4_(iii->fourth(),iii->third(), iii->second(), iii->first()));
      }
    } // end find in cash set
    dih_suite.insert(pair<Spec4, set<Spec4_> >(A,tmpv));
    cout << "." << flush;
  #ifdef CDB
    cout << std::distance(tmpv.begin(),tmpv.end());
  #endif
  } // end every angle
    cout << format("]\n\
  %1$d queries proceeded on database\n") % qry;
    cout << flush;
     */
  }

  /// .. scores is filled with zeroes ..
  void AtomDefiner::scoresZeroFill() {
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;
    qu << format("SELECT `id` FROM `atoms` WHERE `ffield` = %1$d") % atomSettings.ffID;
    res = qu.store();
    // !! removing all scores data !!
    scores.clear();
    assert(res.num_rows());
    map<int, int> zeromap;
    for (co = 0; co < res.num_rows(); ++co) {
      row = res.at(co);
      zeromap.insert(pair<int, int> ( (int) row["id"], 0 ) );
    }
    qu.reset();

    FOR_ATOMS_OF_MOL(it, tp.mol) {
      scores.insert( pair<int, map<int,int> > ( it->GetIdx(), zeromap) );
    }
  } // end scoresZeroFill

  /// .. switching from valence-name to uname ..
  void AtomDefiner::convertAVTtoScores() {
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    map<int, map<int, int> >::iterator score_it;
    map<int, int>::iterator score_subit;
    mysqlpp::Row::size_type co;
    typedef multi_index_container<pair<int, string>,
        indexed_by<
            ordered_unique<
                member<pair<int, string>, int,
                    &pair<int, string>::first> >,
            ordered_non_unique<
                member<pair<int, string>, string,
                    &pair<int, string>::second> > > > unamemap_t;
    // load this array from DB
    unamemap_t mappy;
    qu << format("SELECT id,name FROM atoms WHERE ffield = %1$d") % atomSettings.ffID;
    res = qu.store();
    for (co = 0; co < res.num_rows(); ++co) {
      row = res.at(co);
      mappy.insert(
          pair<int, string>( (int) row["id"], row["name"].c_str() ) );
    }
    qu.reset();

    // apply AVT scores to uname scores using mappy
    for (auto &score: scores) {
      for (auto avtset: avtScores[score.first]) {
        string curname = avtset.first;
        int curscore = avtset.second;
        for (auto iii = mappy.get<1>().lower_bound(curname);
            iii != mappy.get<1>().upper_bound(curname); ++iii) {
          score.second[iii->first] = curscore;
        }
      }
    } // end filling score
  }

  typedef struct {
    int id;
    string type;
    string type2;
    double charge;
    double mass;
    string comment;
  } tempstruct_t;

  void AtomDefiner::atom_align() {
    cout << "Starting atom_alig.." << endl;
    typedef multi_index_container<tempstruct_t,
        indexed_by<
            ordered_unique<member<tempstruct_t, int, &tempstruct_t::id> > > > AtomMapper;
    AtomMapper atom_mapper;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;
    qu
        << "SELECT id, uname, name, charge, mass, comment FROM atoms WHERE ffield = "
        << atomSettings.ffID;
    res = qu.store();
    assert(res);
    cout << "Filling map ..." << flush;
    for (co = 0; co < res.num_rows(); ++co) {
      row = res.at(co);
      tempstruct_t t0;
      t0.id = (int) row["id"];
      t0.type = row["uname"].c_str();
      t0.type2 = row["name"].c_str();
      t0.charge = (double) row["charge"];
      t0.mass = (double) row["mass"];
      t0.comment = row["comment"].c_str();
      atom_mapper.insert(t0);
    }
    assert(co > 1);
    qu.reset();
    cout << "done." << endl;
    // finding atoms with maximum scores
    cout << "Applying scores..." << endl;
    AtomMapper::iterator chk0;
    AtomArray::iterator newa_;
    int max, max_;
    string name;
    for (map<int, map<int, int> >::iterator sit = scores.begin();
        sit != scores.end(); ++sit) {
      max_ = 0;
      for (map<int, int>::iterator mmm = sit->second.begin();
          mmm != sit->second.end(); ++mmm)
        if (mmm->second >= max_) {
          max_ = mmm->second;
          max = mmm->first;
        }
      chk0 = atom_mapper.find(max); // iterator to best atom in atom_mapper
      assert(chk0 != atom_mapper.end());
      name = chk0->type; // name of best atom
      tp.mol.GetAtom(sit->first)->SetType(name);
      AtomArray::iterator newa_ = tp.atoms.find(sit->first);
      assert(newa_ != tp.atoms.end());
      Atom newa = *newa_;
      newa.atom_type = chk0->type;
      newa.atom_type2 = chk0->type2;
      newa.charge = chk0->charge;
      newa.mass = chk0->mass;
      newa.comment = chk0->comment;
      tp.atoms.replace(newa_, newa);
  //    cout << name << endl;
    }
    smart_cgnr();
  }

  /// count scores for znuc, atoms, bonds and dihedrals
  void AtomDefiner::countAVTScores() {
    map<string, int> tmp1, tmp2;
//    map<int, int> tmpit;
    avtScores.clear();
    cout << " ----> Calculating scores for every atom.." << endl;

    // init avtScores, znuc scores apply
    for (auto nb: nbSuite) {
      tmp1.clear();
      for (auto subnb: nb.second)
        tmp1.insert(pair<string, int>(subnb, TPP_ZNUC_COEFF));
      avtScores.insert(pair<int, map<string, int> >(nb.first, tmp1));
    }
    TPPD << "Atomic number map was applied.";

    if (atomSettings.maxbonds) {
      // bond scores apply
      for (auto bond: bondSuite) {
        #ifdef CDB
        cout << "!!!!!!" << bond.first.first() << ":" << bond.first.second() << endl;
        #endif
        for (auto subbond: bond.second) {
          #ifdef CDB
          cout << "===========" << subbond.first() << ":" << subbond.second() << endl;
          #endif
          avtScores[bond.first.first()][subbond.first()] += TPP_BOND_COEFF;
          avtScores[bond.first.second()][subbond.second()] += TPP_BOND_COEFF;
        }
      }
    } // endif maxbond
    TPPD << "Valence bond map was applied.";

    if (atomSettings.maxangles) {
      // angle scores apply
      for (auto angle: angleSuite) {
        #ifdef CDB
        cout << "!!!!!!" << angle.first.first() << ":" << angle.first.second() << ":" << angle.first.third() << endl;
        #endif
        for (auto subangle: angle.second) {
          #ifdef CDB
          cout << "===========" << subangle.first() << ":" << subangle.second() << ":" << subangle.third() << endl;
          #endif
          avtScores[angle.first.first()][subangle.first()] += TPP_ANGLE_COEFF;
          avtScores[angle.first.second()][subangle.second()] += TPP_ANGLE_COEFF;
          avtScores[angle.first.third()][subangle.third()] += TPP_ANGLE_COEFF;
        }
      }
    } // endif angles
    TPPD << "Valence anlge map was applied.";

    if (atomSettings.maxdihedrals) {
      // dihedral scores apply
      for (auto dih: dihdSuite) {
        #ifdef CDB
        cout << "!!!!!!" << dih.first.first() << ":" << dih.first.second() << ":" << dih.first.third() << ":" << dih.first.fourth() << endl;
        #endif
        for (auto subdih: dih.second) {
          #ifdef CDB
          cout << "===========" << subdih.first() << ":" << subdih.second() << ":" << subdih.third() << ":" << subdih.fourth() << endl;
          #endif
          avtScores[dih.first.first()][subdih.first()] += TPP_DIHED_COEFF;
          avtScores[dih.first.second()][subdih.second()] += TPP_DIHED_COEFF;
          avtScores[dih.first.third()][subdih.third()] += TPP_DIHED_COEFF;
          avtScores[dih.first.fourth()][subdih.third()] += TPP_DIHED_COEFF;
        }
      }
    } // end dihedrals
    TPPD << "Dihedral angle map was applied.";

    convertAVTtoScores();
    cout << "                                           ..finished! <---- " << endl;
  }

  void AtomDefiner::printScores(std::ostream &os) {
    // PRINT SCORES
    os << endl;
    for (auto atomscore: scores) {
      os << "FOR ATOM: " << atomscore.first << endl;
      for (auto score: atomscore.second)
        os << format("--- atom # %1$5d: %2$5d pt\n") % score.first
                % score.second;
    }
    os << endl;
  }

  void AtomDefiner::logScores() {
    ostringstream os;
    printScores(os);
    TPPD << os.str();
    cout << "Scores for atom types are written to LOG." << endl;
  }


  /// Standard constructor with DB connection
  AtomDefiner::AtomDefiner(const DbBase::Settings& s1,
                           const AtomDefiner::Settings& s2,
                           Topology &tp_) : DbBase(s1), atomSettings(s2),tp(tp_) {
    connectDB();
    scoresZeroFill();
  }

  /// Standard connection to DB
  bool AtomDefiner::connectDB() {
    DbBase::connectDB();

    ; // what do we need ?

    return true;
  } // connectDB

  /// .. executes sequentially everything that is required for bonded ff ..
  void AtomDefiner::proceed() {
      TPPD << "Filling non-bonded scores.";
      fillNB();
      #ifdef DEBUG
      convertAVTtoScores();
      logScores();
      #endif // DEBUG
      if (atomSettings.maxbonds) {
        TPPD << "Filling bond scores.";
        fillBonds();
        #ifdef DEBUG
        convertAVTtoScores();
        logScores();
        #endif // DEBUG
      }
      if (atomSettings.maxangles) {
        TPPD << "Filling angle scores.";
        fillAngles();
      }
      if (atomSettings.maxdihedrals) {
        TPPD << "Filling dihedral scores.";
        fillDihs();
      }
      convertAVTtoScores();
      smart_fit();
  }


  void AtomDefiner::smart_cgnr() {
        // zero all charge groups
        for (AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
            Atom nat = *it;
            nat.c_gnr = 0;
            tp.atoms.replace(it,nat);
        }
        // algo body
        try {
          TPPD<<"Starting curious SMART-charge-group fitting.";
            mysqlpp::Query qu = con->query();
            qu << format("SELECT id,PAT,flag FROM chargegroups \
                    WHERE ffield = %1$d and flag = 1") % atomSettings.ffID;
            TPPD<<"Loading CGNR patterns from DB..";
            cout << "CHARGEGROUP patterns are loading. Please wait.." << flush;
            QueryResult res;
            res = qu.store();
            if (!res) {
              SqlException e("SQL query failed!");
              e.add("procname", "tpp::AtomDefiner::smart_cgnr");
              e.add("error", "SQL query error");
              e.add("sql_error", qu.error() );
              throw e;
            }
            TPPD<<"OK!\n";
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
                TPPD<<os.str();
                if (sized_patterns.count(pna) == 0)
                    sized_patterns[pna] = vector<string>();
                sized_patterns[pna].push_back(string(row["PAT"]));
                cout << ( format("\b\b\b\b\b%1$4d.") % (int)co ) << flush;
            }
            // Applying patterns from small size to big size
            TppIndex curCG = 1;
            for(std::map<int, vector<string> >::iterator it
                    = sized_patterns.begin(); it != sized_patterns.end();
                    ++it) {
                os.str(""); os.clear();
                os << format("Aplying PATTERNS of size %1$d (%2$d total): \n")
                    % it->first % it->second.size();
                TPPD<<os.str();
                for(vector<string>::iterator ci = it->second.begin();
                        ci != it->second.end(); ++ci) {
                    os.str(""); os.clear();
                    pat.Init(ci->c_str());
                    pat.Match(tp.mol);
                    maplist.clear();
                    maplist = pat.GetUMapList();

                    os << format("[OB] Pattern %1$s matches %2$d times.\n")
                        % (*ci) % maplist.size();
                    TPPD<<os.str();

                    for(int i=0;i<maplist.size();++i) {
                        for (int j=0; j<maplist[i].size(); ++j) {
                            AtomArray::iterator cur_it = tp.atoms.find((int) (maplist[i][j]));
  #ifdef CDB
                            cout << maplist[i][j] << " " << flush;
  #endif
                            assert(cur_it != tp.atoms.end());
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
              std::set<TppIndex> done_atoms;
              std::set<TppIndex> _tempset;
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
                      TppIndex oldcgnr = it->c_gnr;
                      _tempset.clear();
                      for (AtomArray::nth_index_iterator<2>::type cit = tp.atoms.get<2>().lower_bound(oldcgnr);
                          cit != tp.atoms.get<2>().upper_bound(oldcgnr); ++cit) {
                          if (done_atoms.count(cit->index)) continue;
                          _tempset.insert(cit->index);
                          // modifying c_gnr in a separate cycle !! (Index Policy
                          // Needs: see multi_index documentation.. )
                      }
                      for (set<TppIndex>::iterator ii = _tempset.begin(); ii != _tempset.end(); ++ii) {
                          AtomArray::iterator cit = tp.atoms.find(*ii);
                          assert( cit != tp.atoms.end() );
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

        } catch (const Exception &et) {
            cout << "-- CATCH AT SMART_CGNR! --" << endl;
            throw et;
        }
  }

  void AtomDefiner::smart_fit() {

        // make zero-scored copy of scores map
        map<int, map<int, int> > sf_scores (scores);
        for (auto &i : sf_scores)
          for (auto &j : i.second)
            j.second = 0;

        // next work with copied sf_scores
        TPPD<<("Starting curious SMART-fitting procedure.\n");
        mysqlpp::Query qu = con->query();
        QueryResult res;
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
  WHERE  (not atom_patterns.group = 1) and (atoms.ffield = %1$d)") % atomSettings.ffID;
        TPPD<<"Loading patterns from database...";
        cout << "Patterns are loading. Please wait.." << flush;
        res = qu.store();
        if (!res) {
          SqlException e("SQL query failed!");
          e.add("procname", "tpp::AtomDefiner::smartfit");
          e.add("error", "SQL query error");
          e.add("sql_error", qu.error() );
          throw e;
        }
        TPPD<<"OK!";
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
          TPPD<<os.str();
          pat.Match(tp.mol);
          maplist.clear();
          atoms_suite.clear();
          maplist = pat.GetMapList();
  #ifdef CDB
          cout << "============> " << row["PAT"] << "\n"
               << "Matches: " << maplist.size() << endl;
  #endif
          for(int i=0;i<maplist.size();++i) {
            assert( maplist[i].size() >= (int)(row["pos"]) );
            atoms_suite.insert( maplist[i][row["pos"]-1] );
          }

          for(set_it = atoms_suite.begin(); set_it != atoms_suite.end(); ++set_it) {
            score_it = sf_scores.find(*set_it);
            assert( score_it != sf_scores.end() );
            score_subit = score_it->second.find( row["atom_ids"] );
            if ( score_subit == score_it->second.end() ) {
              Exception e("SMARTS-DB ERROR!");
              e.add("procname", "tpp::AtomDefiner::smartfit");
              e.add("error", string("SMART atom pattern #") +
                  lexical_cast<string>(row["apid"]) + " in DB is for invalid atom type!");
              throw e;
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

} // tpp namespace

