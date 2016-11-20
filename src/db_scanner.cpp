#include "mysql.h"
#include "runtime.hpp"
#include "db_scanner.hpp"


#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

//#define CDB

using std::map;
using std::set;
using std::pair;
using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::ostringstream;

using boost::format;
using namespace boost::multi_index;

namespace tpp {
  using namespace OpenBabel;

  bool operator <(const spec2 & a, const spec2 &b) {
    return    (a.first() < b.first()) 
           || ( (a.first() == b.first()) && 
                (a.second() < b.second()) 
              );
  }

  bool operator <(const spec2_ & a, const spec2_ &b) {
    return    (a.first() < b.first()) 
           || ( (a.first() == b.first()) && 
                (a.second() < b.second()) 
              );
  }



  bool operator <(const spec3 & a, const spec3 &b) {
    return (a.first() < b.first()) ||
           ( (a.first() == b.first()) &&
             (spec2(a.second(),a.third()) < spec2(b.second(),b.third())) 
           );           
  }

  bool operator <(const spec3_ & a, const spec3_ &b) {
    return (a.first() < b.first()) ||
           ( (a.first() == b.first()) &&
             (spec2_(a.second(),a.third()) < spec2_(b.second(),b.third())) 
           );           
  }



  bool operator <(const spec4 &a, const spec4 &b) {
    return (a.first() < b.first()) ||
           ( (a.first() == b.first()) &&
             (spec3(a.second(),a.third(),a.fourth()) < spec3(b.second(),b.third(),b.fourth())) 
           );           
  }

  bool operator <(const spec4_ &a, const spec4_ &b) {
    return (a.first() < b.first()) ||
           ( (a.first() == b.first()) &&
             (spec3_(a.second(),a.third(),a.fourth()) < spec3_(b.second(),b.third(),b.fourth())) 
           );           
  }

    void atom_definer::fill_nb() {
      cout << "Comparing atom by elements.." << endl;
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      mysqlpp::Row::size_type co;
      set<int> tmpv;
      map<int, set<int> > znucset;
      int qry = 0;
      cout << "Atoms: [" << flush;
      FOR_ATOMS_OF_MOL(it,tp.mol) {
        int anum = it->GetAtomicNum();
        int ind = it->GetIdx();
#ifdef CDB
        cout << format("#%1$2d,(%2$2d)") % ind % anum;
#endif
        if ( znucset.find(anum) == znucset.end() ) {
          qry++;
          qu = con->query();
          qu << format("\
SELECT MIN(id) as mid FROM atoms \n\
WHERE znuc = %1$d AND ffield = %2$d \n\
GROUP BY name") % anum % ffid;
          res = qu.store();
          BOOST_CHECK(res);
          tmpv.clear();
          for (co = 0; co < res.num_rows(); ++co) {
            row = res.at(co);
            tmpv.insert(row["mid"]);
          }
          znucset.insert(pair<int, set<int> >(anum, tmpv));
        } else {
          tmpv = (znucset.find(anum))->second;
        }
        nb_suite.insert(pair<int,set<int> >(ind, tmpv));
        cout << "." << flush;
        }// atom cycle
      cout << format("]\n\
%1$d queries proceeded on database\n") % qry;
      cout << flush;
      } 

/* ============BONDS====================*/
void atom_definer::fill_bon() {

      cout << "Comparing atom by bonds.." << endl;
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT  res;
      mysqlpp::Row    row;
      mysqlpp::Row::size_type co;
      set<spec2_> tmpv,tmpv1;
      int qry = 0;
      map<spec2_, set<spec2_> > znucset;
      cout << "Bonds: [" << flush;
     FOR_BONDS_OF_MOL(it, tp.mol)
      {
         spec2 A(it->GetBeginAtom()->GetIdx(),it->GetEndAtom()->GetIdx());
         int i0 = A.first(), j0 = A.second();
         int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
             a2 = tp.mol.GetAtom(j0)->GetAtomicNum();
#ifdef CDB
        cout << format("#%1$2d-%2$2d,(%3$2d-%4$2d)") % i0 %j0 % a1 % a2;
#endif
        spec2_ S(a1 <= a2 ? a1 : a2, a1 <= a2 ? a2 : a1);
        if ( znucset.find(S) == znucset.end() ) {
         qry++;
         qu = con->query();
         if (a1 == a2)
         qu << format("\
SELECT \n\
  bonds.id AS bid, \n\
  ( SELECT MIN( id ) FROM atoms WHERE atoms.name = bonds.i AND atoms.ffield = %1$d) AS aid1, \n\
  ( SELECT MIN( id ) FROM atoms WHERE atoms.name = bonds.j AND atoms.ffield = %1$d) AS aid2, \n\
  0 AS inv\n\
FROM bonds\n\
WHERE bonds.ffield = %1$d\n\
AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = bonds.i AND atoms.ffield = %1$d) = %2$d\n\
AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = bonds.j AND atoms.ffield = %1$d) = %2$d\n\
") % this->ffid % a1;
         else
         qu << format("\
( SELECT \n\
  bonds.id AS bid, \n\
  ( SELECT MIN( id ) FROM atoms WHERE atoms.name = bonds.i AND atoms.ffield = %1$d) AS aid1, \n\
  ( SELECT MIN( id ) FROM atoms WHERE atoms.name = bonds.j AND atoms.ffield = %1$d) AS aid2, \n\
  0 AS inv\n\
FROM bonds\n\
WHERE bonds.ffield = %1$d\n\
AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = bonds.i AND atoms.ffield = %1$d) = %2$d\n\
AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = bonds.j AND atoms.ffield = %1$d) = %3$d\n\
) UNION \n\
( SELECT \n\
  bonds.id AS bid, \n\
  ( SELECT MIN( id ) FROM atoms WHERE atoms.name = bonds.i AND atoms.ffield = %1$d) AS aid1, \n\
  ( SELECT MIN( id ) FROM atoms WHERE atoms.name = bonds.j AND atoms.ffield = %1$d) AS aid2, \n\
  1 AS inv\n\
FROM bonds\n\
WHERE bonds.ffield = %1$d\n\
AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = bonds.i AND atoms.ffield = %1$d) = %3$d\n\
AND ( SELECT MIN( znuc ) FROM atoms WHERE atoms.name = bonds.j AND atoms.ffield = %1$d) = %2$d\n\
)\n\
") % this->ffid % a1 % a2;
         res = qu.store();
         BOOST_CHECK(res);
         tmpv.clear();
         for (co = 0; co < res.num_rows(); ++co) {           
            row = res.at(co);
            tmpv.insert((int)row["inv"] == 0 ? 
                    spec2_(row["aid1"],row["aid2"]) : 
                    spec2_(row["aid2"],row["aid1"])
                    );
#ifdef CDB
            cout << "((((((((((" << 
              ((int)row["inv"] == 0 ? row["aid1"] : row["aid2"]) << "---" <<
              ((int)row["inv"] != 0 ? row["aid1"] : row["aid2"]) << endl;

#endif
         }
#ifdef CDB
           cout << "\ninserted: " << S.first() << ":" << S.second() << " elements: " << 
           std::distance(tmpv.begin(),tmpv.end()) <<  endl; 
#endif
         // turning over set
         if (S.first() != a1) { 
           tmpv1.clear();
           for (set<spec2_>::iterator iii = tmpv.begin(); iii != tmpv.end(); ++iii)
             tmpv1.insert(spec2_(iii->second(), iii->first()));
           znucset.insert(pair<spec2_, set<spec2_> >(S, tmpv1));
         } else {
           znucset.insert(pair<spec2_, set<spec2_> >(S, tmpv));
         }
        } else {        
#ifdef CDB
           cout << "FOUND: " << znucset.find(S)->first.first() 
             << "-" << znucset.find(S)->first.second() <<  endl; 
#endif
          tmpv1.clear();
          tmpv1 =  znucset.find(S)->second;
          if (a1 == S.first()) {
           tmpv = tmpv1;
          } else { // turning over set
            tmpv.clear();
            for (set<spec2_>::iterator iii = tmpv1.begin(); iii != tmpv1.end(); ++iii)
             tmpv.insert(spec2_(iii->second(), iii->first()));
          }
        }
        bon_suite.insert(pair<spec2, set<spec2_> >(A,tmpv));
        cout << "." << flush;
#ifdef CDB
        cout << std::distance(tmpv.begin(),tmpv.end());
#endif
      }
      cout << format("]\n\
%1$d queries proceeded on database\n") % qry;
      cout << flush;
    }


  /*=================ANGLES=========================*/
 void atom_definer::fill_ang() {

      cout << "Comparing atom by angles.." << endl;
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      mysqlpp::Row::size_type co;
      set<spec3_> tmpv, tmpv1;
      map<spec3_, set<spec3_> > znucset;
      short qry = 0;
      cout << "Angles: [" << flush;
      FOR_ANGLES_OF_MOL(it, tp.mol)   {
        spec3 A( (*it)[1]+1, (*it)[0]+1, (*it)[2]+1); // one-stabled angle definition
        int i0 = A.first(), j0 = A.second(), k0 = A.third();
        int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
            a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
            a3 = tp.mol.GetAtom(k0)->GetAtomicNum();
#ifdef CDB
       cout << (format("#%1$d-%2$d-%3$d,(%4$d-%5$d-%6$d)") % i0 %j0 % k0 % a1 % a2 % a3) << flush;
#endif
       spec3_ S(a1 <= a3 ? a1 : a3, a2, a1 <= a3 ? a3 : a1);
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
") % this->ffid % a1 % a2;
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
") % this->ffid % a1 % a2 % a3;
        res = qu.store();
        BOOST_CHECK(res);
        tmpv.clear();
        for (co = 0; row = res.at(co); ++co) {
          tmpv.insert((int)row["inv"] == 0 ? 
                       spec3_(row["aid1"],row["aid2"],row["aid3"]) :
                       spec3_(row["aid3"],row["aid2"],row["aid1"])
                     );
        }
#ifdef CDB
           cout << "\ninserted: " << S.first() << ":" << S.second() << ":" << S.third() << " elements: " <<
           std::distance(tmpv.begin(),tmpv.end()) <<  endl; 
#endif

        if (S.first() != a1) { 
           tmpv1.clear();
           for (set<spec3_>::iterator iii = tmpv.begin(); iii != tmpv.end(); ++iii)
             tmpv1.insert(spec3_(iii->third(),iii->second(), iii->first()));
           znucset.insert(pair<spec3_, set<spec3_> >(S, tmpv1));
         } else {
           znucset.insert(pair<spec3_, set<spec3_> >(S, tmpv));
         }
       }  // end find in database
       else {
#ifdef CDB
           cout << "FOUND: " << znucset.find(S)->first.first() 
             << "-" << znucset.find(S)->first.second() 
            << "-" << znucset.find(S)->first.third()  <<  endl; 
#endif
          tmpv1.clear();
          tmpv1 = znucset.find(S)->second;
          if (a1 == S.first()) {
           tmpv = tmpv1;
          } else { // turning over set
            tmpv.clear();
            for (set<spec3_>::iterator iii = tmpv1.begin(); iii != tmpv1.end(); ++iii)
             tmpv.insert(spec3_(iii->third(), iii->second(), iii->first()));
          }
       } // end find in cash set
        ang_suite.insert(pair<spec3, set<spec3_> >(A,tmpv));
        cout << "." << flush;
#ifdef CDB
        cout << std::distance(tmpv.begin(),tmpv.end());
#endif
     } // end every angle
      cout << format("]\n\
%1$d queries proceeded on database\n") % qry;
      cout << flush;
    }



  /*================DIHEDRALS=======================*/
 void atom_definer::fill_dih() {

      cout << "Comparing atom by dihedrals.." << endl;
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      mysqlpp::Row::size_type co;
      set<spec4_> tmpv, tmpv1;
      map<spec4_, set<spec4_> > znucset;
      short qry = 0;
      cout << "Dihedrals: [" << flush;
      FOR_TORSIONS_OF_MOL(it, tp.mol)   {
        spec4 A( (*it)[1]+1, (*it)[0]+1, (*it)[2]+1, (*it)[3]+1); // one-stabled angle definition
        int i0 = A.first(), j0 = A.second(), k0 = A.third(), l0 = A.fourth();
        int a1 = tp.mol.GetAtom(i0)->GetAtomicNum(),
            a2 = tp.mol.GetAtom(j0)->GetAtomicNum(),
            a3 = tp.mol.GetAtom(k0)->GetAtomicNum(),
            a4 = tp.mol.GetAtom(l0)->GetAtomicNum();
#ifdef CDB
       cout << (format("#%1$d-%2$d-%3$d,(%4$d-%5$d-%6$d)") % i0 %j0 % k0 % a1 % a2 % a3) << flush;
#endif
       spec4_ S(a1 <= a4 ? a1 : a4, a1 <= a4 ? a2 : a3, a1 <= a4 ? a3 : a2, a1 <= a4 ? a4 : a1);
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
") % this->ffid % a1 % a2 % a3;
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
") % this->ffid % a1 % a2 % a3 % a4;
//        cout << endl << qu.preview() << endl;
        res = qu.store();
        BOOST_CHECK(res);
        tmpv.clear();
        for (co = 0; row = res.at(co); ++co) {
          tmpv.insert((int)row["inv"] == 0 ? 
                       spec4_(row["aid1"],row["aid2"],row["aid3"],row["aid4"]) :
                       spec4_(row["aid4"],row["aid3"],row["aid2"],row["aid1"])
                     );
        }
#ifdef CDB
           cout << "\ninserted: " << S.first() << ":" << S.second() << ":" << S.third() << " elements: " <<
           std::distance(tmpv.begin(),tmpv.end()) <<  endl; 
#endif

        if (S.first() != a1) { 
           tmpv1.clear();
           for (set<spec4_>::iterator iii = tmpv.begin(); iii != tmpv.end(); ++iii)
             tmpv1.insert(spec4_(iii->fourth(),iii->third(),iii->second(), iii->first()));
           znucset.insert(pair<spec4_, set<spec4_> >(S, tmpv1));
         } else {
           znucset.insert(pair<spec4_, set<spec4_> >(S, tmpv));
         }
       }  // end find in database
       else {
#ifdef CDB
           cout << "FOUND: " << znucset.find(S)->first.first() 
             << "-" << znucset.find(S)->first.second() 
            << "-" << znucset.find(S)->first.third()  <<  endl; 
#endif
          tmpv1.clear();
          tmpv1 = znucset.find(S)->second;
          if (a1 == S.first()) {
           tmpv = tmpv1;
          } else { // turning over set
            tmpv.clear();
            for (set<spec4_>::iterator iii = tmpv1.begin(); iii != tmpv1.end(); ++iii)
             tmpv.insert(spec4_(iii->fourth(),iii->third(), iii->second(), iii->first()));
          }
       } // end find in cash set
        dih_suite.insert(pair<spec4, set<spec4_> >(A,tmpv));
        cout << "." << flush;
#ifdef CDB
        cout << std::distance(tmpv.begin(),tmpv.end());
#endif
     } // end every angle
      cout << format("]\n\
%1$d queries proceeded on database\n") % qry;
      cout << flush;
    }

void atom_definer::spread_atomid() {
      mysqlpp::Query qu = con->query();
      MYSQLPP_RESULT res;
      mysqlpp::Row    row;
      map<int, map<int,int> >::iterator score_it;
      map<int,int>::iterator score_subit;
      mysqlpp::Row::size_type co;
      typedef multi_index_container<
        pair<int,string>,
        indexed_by<
          ordered_unique<member<pair<int,string>, int, &pair<int,string>::first> >,
          ordered_non_unique<member<pair<int,string>, string, &pair<int,string>::second> >
          >
        > unamemap_t;
      unamemap_t mappy;
      qu << format("SELECT id,name FROM atoms WHERE ffield = %1$d") % this->ffid;
      res = qu.store();
      for(co =0; co < res.num_rows(); ++co) {
        row = res.at(co);
        mappy.insert(pair<int,string>((int)row["id"], row["name"].c_str()/*get_string()*/));
      }
      qu.reset();

      for(score_it = scores.begin(); score_it != scores.end(); ++score_it) {
        set<string> nameset;
        map<int,int> tmpmap;
        string curname;
        for(score_subit = score_it->second.begin(); score_subit != score_it->second.end(); 
            ++score_subit) {
          curname = mappy.find(score_subit->first)->second;
          if ( nameset.find(curname) != nameset.end() ) continue;
          nameset.insert(curname);
          for(unamemap_t::nth_index<1>::type::iterator iii = mappy.get<1>().lower_bound(curname);
              iii != mappy.get<1>().upper_bound(curname); ++iii) {
            tmpmap.insert(pair<int,int>(iii->first, score_subit->second));
          }
        }
        score_it->second = tmpmap;
      }
/* SLOW-STYLE spreading...
      for(score_it = scores.begin(); score_it != scores.end(); ++score_it) {
        for(score_subit = score_it->second.begin(); score_subit != score_it->second.end(); ++score_subit) {
          qu << format("\
SELECT atoms.id AS aid FROM atoms \
RIGHT JOIN atoms AS atoms1 ON atoms1.id = %1$d \
WHERE atoms.name = atoms1.name and not (atoms.id = %1$d)") % score_subit->first;
          res = qu.store();
          BOOST_CHECK(res);
          for(co = 0; row = res.at(co); ++co)
            score_it->second.insert(pair<int,int>(row["aid"],score_subit->second));
          qu.reset();
        }
      }
*/
}

  typedef struct { int id; string type; string type2; double charge; double mass; string comment; } tempstruct_t;

void atom_definer::atom_align() {
  cout << "Starting atom_alig.." << endl;
  typedef multi_index_container<
    tempstruct_t,
    indexed_by<
      ordered_unique<member<tempstruct_t, int, &tempstruct_t::id> >
      >
   > AtomMapper;
  AtomMapper atom_mapper;
  mysqlpp::Query qu = con->query();
  MYSQLPP_RESULT res;
  mysqlpp::Row    row;
  mysqlpp::Row::size_type co;
  qu << "SELECT id, uname, name, charge, mass, comment FROM atoms WHERE ffield = " << this->ffid;
  res = qu.store();
  BOOST_CHECK(res);
  cout << "Filling map ..." << flush;
  for (co = 0; co < res.num_rows(); ++co)  {
    row = res.at(co);
    tempstruct_t t0;
    t0.id = (int)row["id"];
    t0.type = row["uname"].c_str();
    t0.type2 = row["name"].c_str();
    t0.charge = (double)row["charge"];
    t0.mass = (double)row["mass"];
    t0.comment = row["comment"].c_str();
    atom_mapper.insert(t0);
  }
  BOOST_CHECK(co > 1);
  qu.reset();
  cout << "done." << endl;
  // finding atoms with maximum scores
  cout << "Applying scores..." << endl;
  AtomMapper::iterator chk0;
  AtomArray::iterator newa_;
  int max, max_; string name;
  for (map<int, map<int,int> >::iterator sit = scores.begin(); sit != scores.end(); ++sit) {
    max_ = 0;
    for (map<int,int>::iterator mmm = sit->second.begin(); mmm != sit->second.end(); ++mmm)
      if (mmm->second >= max_) { max_ = mmm->second; max = mmm->first; }
    chk0 = atom_mapper.find(max); // iterator to best atom in atom_mapper
    BOOST_CHECK(chk0 != atom_mapper.end());
    name = chk0->type; // name of best atom
    tp.mol.GetAtom(sit->first)->SetType(name);
    AtomArray::iterator newa_ = tp.atoms.find(sit->first);
    BOOST_CHECK(newa_ != tp.atoms.end() );
    Atom newa = *newa_;
    newa.atom_type = chk0->type;
    newa.atom_type2 = chk0->type2;
    newa.charge = chk0->charge;
    newa.mass = chk0->mass;
    newa.comment = chk0->comment;
    tp.atoms.replace(newa_,newa);
//    cout << name << endl;
  }
  smart_cgnr();
}

// count scores for znuc, atoms, bonds and dihedrals
void atom_definer::count_scores() {
   map<int,int> tmp, tmp1, tmp2,tmp3;
   map<int,int> tmpit;
   scores.clear();
   cout << "Calculating scores for every atom.." << flush;

   // init scores, znuc scores apply
   for (map<int, set<int> >::iterator it = nb_suite.begin(); it != nb_suite.end(); ++it) {
     tmp.clear();
     for (set<int>::iterator jit=it->second.begin(); jit != it->second.end(); ++jit)
       tmp.insert(pair<int,int>(*jit, TPP_ZNUC_COEFF));
     scores.insert(pair<int, map<int, int> >(it->first, tmp));
   }

   if (PARAM_EXISTS(par, "maxbonds")) {
    // bond scores apply
    for (map<spec2, set<spec2_> >::iterator it = bon_suite.begin(); it != bon_suite.end(); ++it) {
#ifdef CDB
      cout << "!!!!!!" << it->first.first() << ":" << it->first.second() << endl;
#endif
      map<int, map<int,int> >::iterator scit1 = scores.find(it->first.first()),
          scit2 = scores.find(it->first.second());
      BOOST_CHECK( (scit1 != scores.end()) && (scit2 != scores.end()));
      tmp = scit1->second;  // map of 1st atom 
      tmp1 = scit2->second; // map of 2nd atom
      for (set<spec2_>::iterator jit = it->second.begin(); jit != it->second.end(); ++jit) {
#ifdef CDB
        cout << "===========" << jit->first() << ":" << jit->second() << endl;
#endif
        BOOST_CHECK(tmp.find(jit->first()) != tmp.end());
        int old = (tmp.find(jit->first()))->second; // old score
        tmp.erase(jit->first()); // erase from map 
        tmp.insert(pair<int,int>(jit->first(), old+TPP_BOND_COEFF)); // replace new value into map
        BOOST_CHECK(tmp1.find(jit->second()) != tmp1.end());
        old = (tmp1.find(jit->second()))->second;
        tmp1.erase(jit->second());
        tmp1.insert(pair<int,int>(jit->second(), old+TPP_BOND_COEFF));
      }
      scores.erase(it->first.first());
      scores.erase(it->first.second());
      scores.insert(pair<int, map<int,int> >( it->first.first(),  tmp));
      scores.insert(pair<int, map<int,int> >(it->first.second(), tmp1));
    }
   } // endif maxbond

   if (PARAM_EXISTS(par, "maxangles")) {
    // angle scores apply
    for (map<spec3, set<spec3_> >::iterator it = ang_suite.begin(); it != ang_suite.end(); ++it) {
#ifdef CDB
      cout << "!!!!!!" << it->first.first() << ":" << it->first.second() << ":" << it->first.third() << endl;
#endif
      map<int, map<int,int> >::iterator scit1 = scores.find(it->first.first()),
          scit2 = scores.find(it->first.second()),
          scit3 = scores.find(it->first.third());
      BOOST_CHECK( (scit1 != scores.end()) && (scit2 != scores.end()) && (scit3 != scores.end()) );
      tmp = scit1->second;  // map of 1st atom 
      tmp1 = scit2->second; // map of 2nd atom
      tmp2 = scit3->second; // map of 3rd atom
      for (set<spec3_>::iterator jit = it->second.begin(); jit != it->second.end(); ++jit) {
#ifdef CDB
        cout << "===========" << jit->first() << ":" << jit->second() << ":" << jit->third() << endl;
#endif
        BOOST_CHECK(tmp.find(jit->first()) != tmp.end());
        int old = (tmp.find(jit->first()))->second; // old score
        tmp.erase(jit->first()); // erase from map 
        tmp.insert(pair<int,int>(jit->first(), old+TPP_ANGLE_COEFF)); // replace new value into map
        BOOST_CHECK(tmp1.find(jit->second()) != tmp1.end());
        old = (tmp1.find(jit->second()))->second;
        tmp1.erase(jit->second());
        tmp1.insert(pair<int,int>(jit->second(), old+TPP_ANGLE_COEFF));
        BOOST_CHECK(tmp2.find(jit->third()) != tmp2.end());
        old = (tmp2.find(jit->third()))->second;
        tmp2.erase(jit->third());
        tmp2.insert(pair<int,int>(jit->third(), old+TPP_ANGLE_COEFF));
      }
      scores.erase(it->first.first());
      scores.erase(it->first.second());
      scores.erase(it->first.third());
      scores.insert(pair<int, map<int,int> >( it->first.first(),  tmp));
      scores.insert(pair<int, map<int,int> >(it->first.second(), tmp1));
      scores.insert(pair<int, map<int,int> >(it->first.third(), tmp2));
      //insert
    }
   } // endif angles

   if (PARAM_EXISTS(par, "maxdihedrals")) {
    // dihedral scores apply
    for (map<spec4, set<spec4_> >::iterator it = dih_suite.begin(); it != dih_suite.end(); ++it) {
#ifdef CDB
      cout << "!!!!!!" << it->first.first() << ":" << it->first.second() << ":" << it->first.third() << endl;
#endif
      map<int, map<int,int> >::iterator scit1 = scores.find(it->first.first()),
          scit2 = scores.find(it->first.second()),
          scit3 = scores.find(it->first.third()),
          scit4 = scores.find(it->first.fourth());
      BOOST_CHECK(    (scit1 != scores.end()) && (scit2 != scores.end()) 
                    && (scit3 != scores.end()) && (scit4 != scores.end()) );
      tmp = scit1->second;  // map of 1st atom 
      tmp1 = scit2->second; // map of 2nd atom
      tmp2 = scit3->second; // map of 3rd atom
      tmp3 = scit4->second; // map of 3rd atom
      for (set<spec4_>::iterator jit = it->second.begin(); jit != it->second.end(); ++jit) {
#ifdef CDB
        cout << "===========" << jit->first() << ":" << jit->second() << ":" << jit->third() << endl;
#endif
        BOOST_CHECK(tmp.find(jit->first()) != tmp.end());
        int old = (tmp.find(jit->first()))->second; // old score
        tmp.erase(jit->first()); // erase from map 
        tmp.insert(pair<int,int>(jit->first(), old+TPP_DIHED_COEFF)); // replace new value into map
        BOOST_CHECK(tmp1.find(jit->second()) != tmp1.end());
        old = (tmp1.find(jit->second()))->second;
        tmp1.erase(jit->second());
        tmp1.insert(pair<int,int>(jit->second(), old+TPP_DIHED_COEFF));
        BOOST_CHECK(tmp2.find(jit->third()) != tmp2.end());
        old = (tmp2.find(jit->third()))->second;
        tmp2.erase(jit->third());
        tmp2.insert(pair<int,int>(jit->third(), old+TPP_DIHED_COEFF));
        BOOST_CHECK(tmp3.find(jit->fourth()) != tmp3.end());
        old = (tmp3.find(jit->fourth()))->second;
        tmp3.erase(jit->fourth());
        tmp3.insert(pair<int,int>(jit->fourth(), old+TPP_DIHED_COEFF));
      }
      scores.erase(it->first.first());
      scores.erase(it->first.second());
      scores.erase(it->first.third());
      scores.erase(it->first.fourth());
      scores.insert(pair<int, map<int,int> >( it->first.first(),  tmp));
      scores.insert(pair<int, map<int,int> >(it->first.second(), tmp1));
      scores.insert(pair<int, map<int,int> >(it->first.third(), tmp2));
      scores.insert(pair<int, map<int,int> >(it->first.fourth(), tmp3));
      //insert
    }
   // refresh scores on the whole set
   } // endif dihedrals

   spread_atomid();
   cout << " finished!" << endl;
}


void atom_definer::print_scores(std::ostream &os) {
   // PRINT SCORES
   os << endl;
   for (map<int, map<int,int> >::iterator it = scores.begin(); it != scores.end(); ++it) {
     os << "FOR ATOM: " << it->first << endl;
     for (map<int,int>::iterator jit = it->second.begin(); jit != it->second.end(); ++jit)
       os << format("--- atom # %1$5d: %2$5d pt\n") % jit->first % jit->second;
   }
   os << endl;
}

void atom_definer::log_scores() {
  ostringstream os;
  print_scores(os);
  runtime.log_write(os.str());
  cout << "Scores for atom types are written to LOG." << endl;
}

 /*
  * Standard constructor with DB connection
  */
 atom_definer::atom_definer(t_input_params p_, Topology &tp_): db_base(p_), tp(tp_) {
   connect_db();
   this->ffid = boost::lexical_cast<int>(PARAM_READ(this->par, "ffid"));
 }
 
 /*
  * Initial DB queries.
  */
 bool atom_definer::connect_db() {
     db_base::connect_db();

     ; // what do we need ?

     return true;
 } // connect_db

 /*
 * Main class method caller.
 */
 void atom_definer::proceed() {
   try {
     fill_nb();
     if (PARAM_EXISTS(par,"maxbonds"))     fill_bon();
     if (PARAM_EXISTS(par,"maxangles"))    fill_ang();
     if (PARAM_EXISTS(par,"maxdihedrals")) fill_dih();
     count_scores();
     smart_fit();
   } catch(t_sql_exception e) { 
     e.fix_log(); 
     throw e;
   }
 }

} // tpp namespace

