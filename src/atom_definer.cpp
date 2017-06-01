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
    cout << endl;
    TPPI << " ----> Comparing atom by elements..";
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
    cout << "  Atomic numbers: [" << flush;
    for (auto &cur: znucset) {
      qu = con->query();
      ostringstream os;
      os << format(
        "SELECT `atoms`.`name` as nm FROM `atoms` "
        "WHERE `znuc` = %1$d AND `ffield` = %2$d "
        "GROUP BY `name`") % cur.first % atomSettings.ffID;
      #ifdef SQLDEBUG
      TPPD << os.str();
      #endif // SQLDEBUG
      qu << os.str();
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
    mapAItoAVT = znucset;
    cout << "]" << endl;
    TPPD << format("  %1$d queries proceeded on database") % znucset.size();
    #ifdef DEBUG
    ostringstream os;
    os << "\n-- Logging znucset for NB --" << endl;
    for (auto znucelem: znucset) {
      os << znucelem.first << ": ";
      for (auto nm: znucelem.second)
        os << nm << " | ";
      os << endl;
    }
    os << "-- END: logging znucset for NB --" << endl;
    TPPD << os.str();
    #endif // DEBUG
    // fill map to every atom
    FOR_ATOMS_OF_MOL(it, tp.mol) {
      nbSuite.insert(pair<int,set<string> >(it->GetIdx(), znucset[it->GetAtomicNum()]));
    }
    TPPI << "                                   ..finished! <----";
  } // end fillNB

  string AtomDefiner::getSQLSet(int ai) {
    int co = 0;
    string ans = "(";
    for (string e: mapAItoAVT[ai]) {
      if (co) ans += ",";
      ans += "'"+e+"'";
      co++;
    }
    ans += ")";
    return ans;
  }

  /// implement existing bond counting
  void AtomDefiner::fillBonds() {
    TPPI << " ----> Comparing atom by bonds..";
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;

    typedef Spec2_<string> NMBondSpec;   /// bond type specified by atom valence names
    typedef Spec2_<int>    AIBondSpec;   /// bond type specified by atom numbers only
    set<NMBondSpec> tmpv, tmpv1;
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
    cout << "  Bonds: [" << flush;
    for (auto &cur: znucset) {
      int i = cur.first.first(), j = cur.first.second();
        qu = con->query();
        ostringstream os;
        if (i == j)
          os << format(
            " SELECT "
            "   `bonds`.`id` AS `bid`, `bonds`.`i` AS `bi`, `bonds`.`j` AS `bj` "
            " FROM `bonds` "
            " WHERE `bonds`.`ffield` = %1$d "
            "   AND `bonds`.`i` IN %2$s AND `bonds`.`j` IN %2$s "
            ) % atomSettings.ffID % getSQLSet(i);
        else
          os << format(
            " SELECT "
            "       `bonds`.`id` AS `bid`, `bonds`.`i` AS `bi`, `bonds`.`j` AS `bj` "
            "     FROM bonds "
            "     WHERE `bonds`.`ffield` = %1$d "
            "       AND `bonds`.`i` IN %2$s AND `bonds`.`j` IN %3$s "
            "     "
            " UNION "
            "     "
            " SELECT "
            "       `bonds`.`id` AS `bid`, `bonds`.`j` AS `bi`, `bonds`.`i` AS `bj` "
            "     FROM bonds "
            "     WHERE `bonds`.`ffield` = %1$d "
            "       AND `bonds`.`i` IN %3$s AND `bonds`.`j` IN %2$s "
            ) % atomSettings.ffID % getSQLSet(i) % getSQLSet(j);
        #ifdef SQLDEBUG
        TPPD << os.str();
        #endif // SQLDEBUG
        qu << os.str();
        res = qu.store();
        if (!res) {
          Exception e("Empty bond list resulted.");
          e.add("query", qu.str());
          throw e;
        }
        for (co = 0; co < res.num_rows(); ++co) {
          row = res.at(co);
          cur.second.insert( NMBondSpec(string(row["bi"]),string(row["bj"]) ) );
        }
        cout << "." << flush;
    } // end for znucset
    cout << "]" << endl;
    TPPI << format("  %1$d queries proceeded on database") % znucset.size();
    #ifdef DEBUG
    ostringstream os;
    TPPD << "-- Logging znucset for BONDS --";
    for (auto znucelem: znucset) {
      os.str("");
      os << format("(%1$d-%2$d): ") % znucelem.first.first() % znucelem.first.second();
      for (auto nmpair: znucelem.second) {
        os << format("%1$s-%2$s | ") % nmpair.first() % nmpair.second();
      }
      TPPD << os.str();
    }
    TPPD << "-- END: Logging znucset for BONDS --";
    #endif // DEBUG

    // associate znucset bondsets with every bond
    FOR_BONDS_OF_MOL(it, tp.mol) {
      Spec2<int> idxBnd(it->GetBeginAtom()->GetIdx(), it->GetEndAtom()->GetIdx()); // bond i->j
      int i0 = idxBnd.first(), j0 = idxBnd.second(); // i,j
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum();
      TPPD << format("Inserting %d - %d") % idxBnd.first() % idxBnd.second();
      if ( znucset.count(AIBondSpec(a1,a2)) ) {
        bondSuite.insert( pair<Spec2<int>, set<NMBondSpec> > (
           idxBnd, znucset[AIBondSpec(a1,a2)] ) );
      } else if ( znucset.count(AIBondSpec(a2,a1) ) ) {
        // swap first and second in the set
        tmpv1.clear();
        for (auto ii: znucset[AIBondSpec(a2,a1)])
          tmpv1.insert(NMBondSpec(ii.second(),ii.first()));
        bondSuite.insert( pair<Spec2<int>, set<NMBondSpec> > (
           idxBnd, set<NMBondSpec>(tmpv1) ) );
      } else {
        /// This situation can not occur. Just test that the code is correct.
        TPPE << "No AIBondSpec was prepared!";
        assert(0);
      }
    } // end FOR_BONDS_OF_MOL
    TPPI << "                                ..finished! <---- ";
  } // end fillBonds

  /// implements angle counting
  void AtomDefiner::fillAngles() {
    cout << "Comparing atom by angles.." << endl;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;

    typedef Spec3_<string> NMAngleSpec;   /// angle type specified by atom valence names
    typedef Spec3_<int>    AIAngleSpec;   /// angle type specified by atom numbers only
    set<NMAngleSpec> tmpv, tmpv1;
    int qry = 0;
    map<AIAngleSpec, set<NMAngleSpec> > znucset; /// map: AtomNumber-AtomNumber-AtomNumber -> AtomName-AtomName-AtomName

    // all angle types are included once
    FOR_ANGLES_OF_MOL(it, tp.mol) {
      int i0 = (*it)[1] + 1, j0 = (*it)[0] + 1, k0 = (*it)[2] + 1; // i,j,k
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
          a3 = tp.mol.GetAtom(k0)->GetAtomicNum();
      AIAngleSpec S(a1 <= a3 ? a1 : a3, a2, a1 <= a3 ? a3 : a1);
      znucset.insert( pair<AIAngleSpec, set<NMAngleSpec> > ( S, set<NMAngleSpec>() ) );
    }

    cout << "Angles: [" << flush;
    for (auto &cur: znucset) {
      int i = cur.first.first(), j = cur.first.second(), k = cur.first.third();
        qu = con->query();
        ostringstream os;
        if (i == k)
          os << format(
            " SELECT "
            "   `angles`.`id` AS `aid`, `angles`.`i` AS `ai`, `angles`.`j` AS `aj`, `angles`.`k` AS `ak` "
            " FROM `angles` "
            " WHERE `angles`.`ffield` = %1$d "
            "   AND `angles`.`i` IN %2$s AND `angles`.`j` IN %3$s AND `angles`.`k` IN %2$s "
            ) % atomSettings.ffID % getSQLSet(i) % getSQLSet(j);
        else
          os << format(
            " SELECT "
            "     `angles`.`id` AS `aid`, `angles`.`i` AS `ai`, `angles`.`j` AS `aj`, `angles`.`k` AS `ak` "
            "   FROM `angles` "
            "   WHERE `angles`.`ffield` = %1$d "
            "     AND `angles`.`i` IN %2$s AND `angles`.`j` IN %3$s AND `angles`.`k` IN  %4$s "
            "     "
            " UNION "
            "     "
            " SELECT "
            "     `angles`.`id` AS `aid`, `angles`.`k` AS `ai`, `angles`.`j` AS `aj`, `angles`.`i` AS `ak` "
            "   FROM `angles` "
            "   WHERE `angles`.`ffield` = %1$d "
            "     AND `angles`.`i` IN %4$s AND `angles`.`j` IN %3$s AND `angles`.`k` IN  %2$s "
            ) % atomSettings.ffID % getSQLSet(i) % getSQLSet(j) % getSQLSet(k);
        #ifdef SQLDEBUG
        TPPD << os.str();
        #endif // SQLDEBUG
        qu << os.str();
        res = qu.store();
        if (!res) {
          Exception e("Empty angle list resulted.");
          e.add("query", qu.str());
          throw e;
        }
        for (co = 0; co < res.num_rows(); ++co) {
          row = res.at(co);
          cur.second.insert( NMAngleSpec(string(row["ai"]),string(row["aj"]), string(row["ak"]) ) );
        }
        cout << "." << flush;
    } // end for znucset
    cout << format("]\n"
              "%1$d queries proceeded on database\n") % znucset.size();
    cout << flush;
    TPPD << format("%1$d queries proceeded on database\n") % znucset.size();

    #ifdef DEBUG
    ostringstream os;
    os << "\n-- Logging znucset for ANGLES --" << endl;
    for (auto znucelem: znucset) {
      os << format("(%1$d-%2$d-%3$d): ") % znucelem.first.first() % znucelem.first.second() % znucelem.first.third();
      for (auto nmtrt: znucelem.second) {
        os << format("%1$s-%2$s-%3$s | ") % nmtrt.first() % nmtrt.second() % nmtrt.third();
      }
      os << endl;
    }
    os << "\n-- END: Logging znucset for ANGLES --" << endl;
    TPPD << os.str();
    #endif // DEBUG

    // associate znucset bondsets with every angle
    FOR_ANGLES_OF_MOL(it, tp.mol) {
      Spec3<int> idxAng( (*it)[1] + 1, (*it)[0] + 1, (*it)[2] + 1); // angle i->j<-k
      int i0 = idxAng.first(), j0 = idxAng.second(), k0 = idxAng.third(); // i,j,k
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
          a3 = tp.mol.GetAtom(k0)->GetAtomicNum();
      if ( znucset.count(AIAngleSpec(a1,a2,a3)) ) {
        angleSuite.insert( pair<Spec3<int>, set<NMAngleSpec> > (
           idxAng, znucset[AIAngleSpec(a1,a2,a3)] ) );
      } else if ( znucset.count(AIAngleSpec(a3,a2,a1) ) ) {
        // swap first and second in the set
        tmpv1.clear();
        for (auto ii: znucset[AIAngleSpec(a3,a2,a1)])
          tmpv1.insert(NMAngleSpec(ii.third(),ii.second(),ii.first()));
        angleSuite.insert( pair<Spec3<int>, set<NMAngleSpec> > (
           idxAng, tmpv1 ) );
      } else {
        /// This situation can not occur. Just test that the code is correct.
        TPPE << format("No AIAngleSpec was prepared: %d-%d-%d ") % a1 % a2 % a3 ;
        assert(0);
      }
    } // end FOR_ANGLES_OF_MOL

  } // end fillAngles

  /// implements dihedral counting
  void AtomDefiner::fillDihs() {
    cout << "Comparing atom by dihedrals.." << endl;
    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;

    typedef Spec4_<string> NMDihdSpec;   /// angle type specified by atom valence names
    typedef Spec4_<int>    AIDihdSpec;   /// angle type specified by atom numbers only
    set<NMDihdSpec> tmpv, tmpv1;
    int qry = 0;
    map<AIDihdSpec, set<NMDihdSpec> > znucset; /// map: AtomNumber-AN-AN-AN -> AtomName-ANm-ANm-ANm

    // all angle types are included once
    FOR_TORSIONS_OF_MOL(it, tp.mol) {
      int i0 = (*it)[0] + 1, j0 = (*it)[1] + 1, k0 = (*it)[2] + 1, l0 = (*it)[3] + 1; // i,j,k
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
          a3 = tp.mol.GetAtom(k0)->GetAtomicNum(),
          a4 = tp.mol.GetAtom(l0)->GetAtomicNum();
      bool ori = ( (a1 == a4) && (a2 > a3) ) || (a1 > a4);
      AIDihdSpec S(ori ? a1 : a4, ori ? a2 : a3, ori ? a3 : a2, ori ? a4 : a1);
      znucset.insert( pair<AIDihdSpec, set<NMDihdSpec> > ( S, set<NMDihdSpec>() ) );
    }

    cout << "Dihedrals: [" << flush;
    for (auto &cur: znucset) {
      int i = cur.first.first(), j = cur.first.second(), k = cur.first.third(), l = cur.first.fourth();
        qu = con->query();
        ostringstream os;
        if ( (i == l) and (j == k) )
          os << format(
            " SELECT "
            "   `dihedrals`.`id` AS `did`, `dihedrals`.`i` AS `di`, `dihedrals`.`j` AS `dj`, `dihedrals`.`k` AS `dk`, `dihedrals`.`l` AS `dl` "
            " FROM `dihedrals` "
            " WHERE `dihedrals`.`ffield` = %1$d "
            "   AND `dihedrals`.`i` IN %2$s AND `dihedrals`.`j` IN %3$s AND `dihedrals`.`k` IN %3$s AND `dihedrals`.`l` IN %2$s "
            ) % atomSettings.ffID % getSQLSet(i) % getSQLSet(j);
        else
          os << format(
            " SELECT "
            "     `dihedrals`.`id` AS `did`, `dihedrals`.`i` AS `di`, `dihedrals`.`j` AS `dj`, `dihedrals`.`k` AS `dk`, `dihedrals`.`l` AS `dl` "
            "   FROM `dihedrals` "
            "   WHERE `dihedrals`.`ffield` = %1$d "
            "     AND `dihedrals`.`i` IN %2$s AND `dihedrals`.`j` IN %3$s AND `dihedrals`.`k` IN %4$s AND `dihedrals`.`l` IN %5$s "
            "     "
            " UNION "
            "     "
            " SELECT "
            "     `dihedrals`.`id` AS `did`, `dihedrals`.`l` AS `di`, `dihedrals`.`k` AS `dj`,  `dihedrals`.`j` AS `dk`, `dihedrals`.`i` AS `dl`"
            "   FROM `dihedrals` "
            "   WHERE `dihedrals`.`ffield` = %1$d "
            "     AND `dihedrals`.`i` IN %5$s AND `dihedrals`.`j` IN %4$s AND `dihedrals`.`k` IN %3$s AND `dihedrals`.`l` IN %2$s "
            ) % atomSettings.ffID % getSQLSet(i) % getSQLSet(j) % getSQLSet(k) % getSQLSet(l);
        #ifdef SQLDEBUG
        TPPD << os.str();
        #endif // SQLDEBUG
        qu << os.str();
        res = qu.store();
        if (!res) {
          Exception e("Empty dihedral list resulted.");
          e.add("query", qu.str());
          throw e;
        }
        for (co = 0; co < res.num_rows(); ++co) {
          row = res.at(co);
          cur.second.insert( NMDihdSpec(string(row["di"]),string(row["dj"]), string(row["dk"]), string(row["dl"]) ) );
        }
        cout << "." << flush;
    } // end for znucset
    cout << format("]\n"
              "%1$d queries proceeded on database\n") % znucset.size();
    cout << flush;
    TPPD << format("%1$d queries proceeded on database\n") % znucset.size();

    #ifdef DEBUG
    ostringstream os;
    os << "\n-- Logging znucset for DIHEDRALS --" << endl;
    for (auto znucelem: znucset) {
      os << format("(%1$d-%2$d-%3$d-%4$d): ") % znucelem.first.first() % znucelem.first.second() % znucelem.first.third() % znucelem.first.fourth();
      for (auto nmqd: znucelem.second) {
        os << format("%1$s-%2$s-%3$s-%4$s | ") % nmqd.first() % nmqd.second() % nmqd.third() % nmqd.fourth();
      }
      os << endl;
    }
    os << "\n-- END: Logging znucset for DIHEDRALS --" << endl;
    TPPD << os.str();
    #endif // DEBUG

    // associate znucset bondsets with every angle
    FOR_TORSIONS_OF_MOL(it, tp.mol) {
      Spec4<int> idxDih( (*it)[0] + 1, (*it)[1] + 1, (*it)[2] + 1, (*it)[3] + 1); // angle i-j-k-l
      int i0 = idxDih.first(), j0 = idxDih.second(), k0 = idxDih.third(), l0 = idxDih.fourth(); // i,j,k
      int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
          a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
          a3 = tp.mol.GetAtom(k0)->GetAtomicNum(),
          a4 = tp.mol.GetAtom(l0)->GetAtomicNum();
      if ( znucset.count(AIDihdSpec(a1,a2,a3,a4)) ) {
        dihdSuite.insert( pair<Spec4<int>, set<NMDihdSpec> > (
           idxDih, znucset[AIDihdSpec(a1,a2,a3,a4)] ) );
      } else if ( znucset.count(AIDihdSpec(a4,a3,a2,a1) ) ) {
        // swap first and second in the set
        tmpv1.clear();
        for (auto ii: znucset[AIDihdSpec(a4,a3,a2,a1)])
          tmpv1.insert(NMDihdSpec(ii.fourth(),ii.third(),ii.second(),ii.first()));
        dihdSuite.insert( pair<Spec4<int>, set<NMDihdSpec> > (
           idxDih, tmpv1 ) );
      } else {
        /// This situation can not occur. Just test that the code is correct.
        TPPE << format("No AIDihdSpec was prepared: %d-%d-%d-%d ") % a1 % a2 % a3 % a4 ;
        assert(0);
      }
    } // end FOR_TORSIONS_OF_MOL

  }

  /// .. scores is filled with zeroes ..
  void AtomDefiner::scoresZeroFill() {
    scores.clear();   // !! removing all scores data !!

    map<int, int> zeromap;
    for (auto it = atom_mapper.begin(); it != atom_mapper.end(); ++it) {
      zeromap.insert(pair<int, int> ( it->id, 0 ) );
    }

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

  /// general atomtype attribution
  void AtomDefiner::atomAlign() {
    TPPI << " ----> Finalization of atomtype attribution..";
    // finding atoms with maximum scores
    TPPI << "  Applying scores.";
    AtomMapper::iterator chk0;
    AtomArray::iterator newa_;
    int max, max_;
    string name;
    for (auto sit: scores) {
      max_ = 0;
      for (auto mmm: sit.second) /// find type with max scores
        if (mmm.second >= max_) {
          max_ = mmm.second;
          max = mmm.first;
        }
      chk0 = atom_mapper.find(max); // iterator to best atom in atom_mapper
      assert(chk0 != atom_mapper.end());
      name = chk0->type; // name of best atom
      tp.mol.GetAtom(sit.first)->SetType(name);
      AtomArray::iterator newa_ = tp.atoms.find(sit.first);
      assert(newa_ != tp.atoms.end());
      Atom newa = *newa_;
      newa.atom_type = chk0->type;
      newa.atom_type2 = chk0->type2;
      newa.charge = chk0->charge;
      newa.mass = chk0->mass;
      newa.comment = chk0->comment;
      tp.atoms.replace(newa_, newa);
    }

    smartCgnr();

    TPPI << "                                             ..DONE. <----";
  } // end atomAlign procedure

  /// count scores for znuc, atoms, bonds and dihedrals
  void AtomDefiner::countAVTScores() {
    map<string, int> tmp1, tmp2;
    TPPD << "Valence-type scores (avtScores) are cleared.";
    avtScores.clear();
    cout << endl;
    TPPI << " ----> Calculating scores for every atom..";

    // init avtScores, znuc scores apply
    for (auto nb: nbSuite) {
      tmp1.clear();
      for (auto subnb: nb.second)
        tmp1.insert(pair<string, int>(subnb, TPP_ZNUC_COEFF));
      avtScores.insert(pair<int, map<string, int> >(nb.first, tmp1));
    }
    TPPD << "  Atomic number map was applied.";

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
    TPPD << "  Valence bond map was applied.";

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
    TPPD << "  Valence anlge map was applied.";

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
    TPPD << "  Dihedral angle map was applied.";

    convertAVTtoScores();
    TPPI << "                                           ..finished! <---- ";
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
    cout << endl;
    TPPI << "== Starting AtomDefiner ==";
    checkAtomlistConsistency();
    connectDB();
    scoresZeroFill();
  }

  /// crucial internal check
  void AtomDefiner::checkAtomlistConsistency() {
    FOR_ATOMS_OF_MOL(it,tp.mol) {
      AtomArray::iterator pa = tp.atoms.find(it->GetIdx());

      // self-consistency check
      if ( pa == tp.atoms.end() ) {
        Exception e("OBMol - tpp::AtomArray inconsistensy!");
        e.add("procname", "tpp::AtomDefiner::checkAtomlistConsistency");
        e.add("error", "There are atom which is absent in tp.atoms: "+boost::lexical_cast<string>(it->GetIdx()) );
        throw e;
      }
    }
  }

  /// Standard connection to DB
  bool AtomDefiner::connectDB() {
    DbBase::connectDB();

    mysqlpp::Query qu = con->query();
    QueryResult res;
    mysqlpp::Row row;
    mysqlpp::Row::size_type co;
    ostringstream os;

    TPPD << "  Requesting full atomtype table.";
    os << "SELECT id, uname, name, charge, mass, comment FROM atoms WHERE ffield = "
        << atomSettings.ffID;
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    qu << os.str();
    res = qu.store();
    if (!res) {
      SqlException e("Error in SQL request.");
      e.add("procname", "tpp::AtomDefiner::atomAlign");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error() );
      e.add("query", os.str());
      throw e;
    }
    TPPD << "  Forming atomtype map.";
    // add UNDEF atom
    {
      tempstruct_t t0;
      t0.id = -1;
      t0.type = "undef";
      t0.type2 = "UU";
      t0.charge = 0.0;
      t0.mass = 0.0;
      t0.comment = "Undefined type [has no sense]";
      atom_mapper.insert(t0);
    }
    // add normat atomtype
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
    if (!co) {
      SqlException e("Empty atomtype list resulted.");
      e.add("procname", "tpp::AtomDefiner::atomAlign");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error() );
      e.add("query", os.str());
      throw e;
    }
    qu.reset();

    return true;
  } // connectDB

  /// .. executes sequentially everything that is required for bonded ff ..
  void AtomDefiner::proceed() {
      TPPD << "Filling non-bonded scores.";
      fillNB();
      #ifdef DEBUG
      countAVTScores();
      logScores();
      #endif // DEBUG
      if (atomSettings.maxbonds) {
        TPPD << "Filling bond scores.";
        fillBonds();
        #ifdef DEBUG
        countAVTScores();
        logScores();
        #endif // DEBUG
      }
      if (atomSettings.maxangles) {
        TPPD << "Filling angle scores.";
        fillAngles();
        #ifdef DEBUG
        countAVTScores();
        logScores();
        #endif // DEBUG
      }
      if (atomSettings.maxdihedrals) {
        TPPD << "Filling dihedral scores.";
        fillDihs();
        #ifdef DEBUG
        countAVTScores();
        logScores();
        #endif // DEBUG
      }
      countAVTScores();
      smartFit();
      #ifdef DEBUG
      logScores();
      #endif // DEBUG
  }

  /// implementing charge group definition
  void AtomDefiner::smartCgnr() {
        ostringstream os;
        // zero all charge groups
        for (AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
            Atom nat = *it;
            nat.c_gnr = 0;
            tp.atoms.replace(it,nat);
        }
        // algo body
        try {
            TPPI << "  Starting charge-group attribution algorithm.";
            mysqlpp::Query qu = con->query();
            os << format("SELECT id,PAT,flag FROM chargegroups "
                    "WHERE ffield = %1$d and flag = 1") % atomSettings.ffID;
            #ifdef SQLDEBUG
            TPPD << os.str();
            #endif
            qu << os.str();
            TPPD << "  Loading CGNR patterns from DB.";
            QueryResult res;
            res = qu.store();
            if (!res) {
              SqlException e("SQL query failed!");
              e.add("procname", "tpp::AtomDefiner::smartCgnr");
              e.add("error", "SQL query error");
              e.add("sql_error", qu.error() );
              e.add("query", os.str());
              throw e;
            }

            OBSmartsPattern pat;
            mysqlpp::Row    row;
            mysqlpp::Row::size_type co;
            vector<vector<int> > maplist;
            set<int> atoms_suite;
            std::map<int, vector<string> > sized_patterns;
            std::ostringstream os;
            TPPD << " Matching CGNR patterns.";
            if (!atomSettings.verbose)
              cout << ( format("  Patterns checked: %1$4d.") % 0 ) << flush;
            for (co=0; co < res.num_rows(); ++co) {
                row = res.at(co);
                pat.Init(row["PAT"]);
                int pna = pat.NumAtoms();
                TPPD << format("  ** [OB] Process PAT: %1$s having %2$d atoms.")
                    % row["PAT"] % pna;
                if (sized_patterns.count(pna) == 0)
                    sized_patterns[pna] = vector<string>();
                sized_patterns[pna].push_back(string(row["PAT"]));
                if (!atomSettings.verbose)
                  cout << ( format("\b\b\b\b%1$4d") % (int)co ) << flush;
            }

            if (!atomSettings.verbose)
              cout << ( format("\n  Patterns applied: %1$4d.") % 0 ) << flush;
            // Applying patterns from small size to big size
            TppIndex curCG = 1;
            for(std::map<int, vector<string> >::iterator it
                    = sized_patterns.begin(); it != sized_patterns.end();
                    ++it) {
                TPPD << format("  Aplying PATTERNS of size %1$d (%2$d total):")
                    % it->first % it->second.size();
                for(vector<string>::iterator ci = it->second.begin();
                        ci != it->second.end(); ++ci) {
                    pat.Init(ci->c_str());
                    pat.Match(tp.mol);
                    maplist.clear();
                    maplist = pat.GetUMapList();

                    TPPD << format("  ** [OB] Pattern %1$s matches %2$d times.")
                        % (*ci) % maplist.size();

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
                    if (!atomSettings.verbose)
                      cout << ( format("\b\b\b\b%1$4d") % (int)(curCG-1) ) << flush;
                } // cycle of paterns of equal size
            } //cycle over sizes
            if (!atomSettings.verbose)
              cout << endl;

            { // independent block of CGR renumbering
              TPPI << "  Renumbering CGNR according to human-readable style.";
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
            } // ending CGNR renumbering block

        } catch (const Exception &et) {
            cout << "-- CATCH AT SMART_CGNR! --" << endl;
            throw et;
        }
  }

  /// implementing fit according to SMARTS db
  void AtomDefiner::smartFit() {
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
        ostringstream os;

        // make zero-scored copy of scores map
        map<int, map<int, int> > sf_scores (scores);
        for (auto &i : sf_scores)
          for (auto &j : i.second)
            j.second = 0;
        // make map of fitted patterns & fill it with empty strings
        map<int, map<int, string> > sf_scores_smarts;
        for (auto &i: sf_scores) {
          map<int,string> m0;
          for (auto &j: i.second) {
            m0.insert(pair<int,string>(j.first, ""));
          }
          sf_scores_smarts.insert(pair<int, map<int,string> >(i.first, m0) );
        }


        // next work with copied sf_scores
        cout << endl;
        TPPI << " ----> Performing SMART-based atomtype attribution..";
        os << format("\
  SELECT atom_patterns.id as apid, PAT, pos, atom_ids, atoms.znuc AS znuc, good \
  FROM atom_patterns \
  RIGHT JOIN atoms ON atoms.id = atom_patterns.atom_ids \
  WHERE  (not atom_patterns.group = 1) and (atoms.ffield = %1$d)") % atomSettings.ffID;
        #ifdef SQLDEBUG
        TPPD << os.str();
        #endif // SQLDEBUG
        qu << os.str();
        TPPI << "  Loading patterns from database.";
        res = qu.store();
        if (!res) {
          SqlException e("SQL query failed!");
          e.add("procname", "tpp::AtomDefiner::smartfit");
          e.add("error", "SQL query error");
          e.add("sql_error", qu.error() );
          throw e;
        }
        TPPD << "  SMART patterns have been loaded from DB.";
        TPPI << "  Checking every SMART-pattern.";
        if (!atomSettings.verbose)
          cout << ( format("  Patterns checked: %1$4d.") % 0 ) << flush;

        // process every pattern while reading DB rows
        for (co=0; co < res.num_rows(); ++co) {
          row = res.at(co);
          pat.Init(row["PAT"]);
          std::ostringstream os;
          os << format("  ** [OB] Process PAT: %1$s having %2$d atoms.")
              % row["PAT"] % pat.NumAtoms();
          TPPD << os.str();
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

          for (set_it = atoms_suite.begin(); set_it != atoms_suite.end(); ++set_it) {
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
              os.str("");
              os << format(" <=== Good: %-4d No %-2d Pat \"%s\"") % (int)row["good"] % row["pos"] % row["PAT"] << flush;
              sf_scores_smarts[*set_it][row["atom_ids"]] = os.str();
            }

          }
          if (!atomSettings.verbose)
            cout << ( format("\b\b\b\b\b%1$4d.") % (int)co ) << flush;
        } // scanning over every smarts

        if (!atomSettings.verbose)
          cout << "\n";

        printSmartFitStats(sf_scores, sf_scores_smarts);

        // check for SMART was found for every atom
        int bad_atom_num = 0;
        for (auto &i : sf_scores) {
          auto maxptr = std::max_element(
            std::begin(i.second), std::end(i.second),
            [] (const pair<int,int> & p1, const pair<int,int> & p2) {
                return p1.second < p2.second;
            }
          );
          if (maxptr->second < 50) { // TODO: define 50 as settable threshold
            bad_atom_num++;
          }
        }
        if (bad_atom_num) {
          TPPI << format("Atoms not normally recognized with SMARTS: %d") % bad_atom_num;

          if (!atomSettings.maxbonds && !atomSettings.maxangles && !atomSettings.maxdihedrals) {
              Exception e("UNRECOGNIZED ATOMS IN YOUR SYSTEM!");
              e.add("procname", "tpp::AtomDefiner::smartfit");
              e.add("error",
                "Your molecule have atoms in unrecognized chemical environment. "
                "See SMART FIT STATISTICS in your log file for details. "
                "You can not proceed topology construction without additional atomtype searching steps. "
                "It is recomended to add unrecognized type to DB. "
                "If you are not going to modify DB, please try to rerun with --max-bonds. ");
              throw e;
          }
        }

        // apply smart scores to main scores map
        for (auto &i : sf_scores)
          for (auto &j : i.second)
            scores[i.first][j.first] += TPP_SMART_COEFF * j.second;

        // I << " ----> Performing SMART-based atomtype attribution..";
        TPPI << "                                                    ..DONE. <----";
  } // end smartfit

  /// print SMART fit statistics
  void AtomDefiner::printSmartFitStats(map<int, map<int,int> > &_sfscores,
    map<int,map<int,string> > &_sfsmarts) {
      TPPD << "Generation of SMART statistics started.\n\n";
      TPPD << "                   SMART FIT STATISTICS\n";
      TPPD << "========================================================================";
      TPPD << " id | name  |  uname        | v.name  | smart score|pos in pt.| pattern ";
      TPPD << "------------------------------------------------------------------------";
      // the following strings are generated in report:
             //   1: C1   ==> opls_145        CA    <=== Good: 140  No 1  Pat "[cr6]"
      // finding the best one!
      for (auto i : _sfscores) {
        auto maxptr = std::max_element(
          std::begin(i.second), std::end(i.second),
          [] (const pair<int,int> & p1, const pair<int,int> & p2) {
              return p1.second < p2.second;
          }
        );
        string smartStr = _sfsmarts[i.first][maxptr->first];
        auto atit = atom_mapper.find(maxptr->first);
        auto atom = tp.atoms.find(i.first);
        TPPD << format("%4d: %-4s ==> %-15s %-4s %s")
         % i.first % atom->atom_name % atit->type % atit->type2 % smartStr;
      }
      TPPD << "========================================================================\n\n";
  } // end printSmartFitStats

} // tpp namespace

