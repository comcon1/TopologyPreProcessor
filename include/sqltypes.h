/*
 * ATTENTION! 
 *
 * Because of SSQL macros you can not include this header into another header!
 * Include it to .cpp directly!
 */

#ifndef TPP_SQLTYPES_H
#define TPP_SQLTYPES_H

#include <mysql++/mysql++.h>
#include <mysql++/ssqls.h>

  sql_create_5(selftest,
      1,5,
      mysqlpp::sql_int, id,
      mysqlpp::sql_int, ffid,
      mysqlpp::sql_text, pdb,
      mysqlpp::sql_int, charge,
      mysqlpp::sql_varchar, comment)

  sql_create_15(selftestcase,
      1,15,
      mysqlpp::sql_int,        id,
      mysqlpp::sql_int,        molid,
      mysqlpp::sql_int,        ffid,
      mysqlpp::sql_enum,       type1,
      mysqlpp::sql_varchar,    type2,
      mysqlpp::sql_int,        i,
      mysqlpp::sql_int_null,   j,
      mysqlpp::sql_int_null,   k,
      mysqlpp::sql_int_null,   l,
      mysqlpp::sql_float,      c1,
      mysqlpp::sql_float_null, c2,
      mysqlpp::sql_float_null, c3,
      mysqlpp::sql_float_null, c4,
      mysqlpp::sql_float_null, c5,
      mysqlpp::sql_float_null, c6)


namespace tpp {

  class selftestcase_Atom : public selftestcase {
    public:
      selftestcase_Atom(const mysqlpp::sql_int &molid, const mysqlpp::sql_int &ffid, const t_atom& a) :
        selftestcase( 0, molid, ffid, 
            "atom", a.atom_type, a.index, mysqlpp::null, mysqlpp::null, mysqlpp::null,
            a.charge, a.mass, mysqlpp::null, mysqlpp::null, mysqlpp::null, mysqlpp::null ) {

      }
  };


  class selftestcase_Bond : public selftestcase {
    public:
      selftestcase_Bond(const mysqlpp::sql_int &molid, const mysqlpp::sql_int &ffid, const TopCoord& c, const TopElement &e) :
        selftestcase( 0, molid, ffid, 
            "bond", boost::lexical_cast<string>(c.f), e.i, e.j, mysqlpp::null, mysqlpp::null,
            c.c0, c.c1, mysqlpp::null, mysqlpp::null, mysqlpp::null, mysqlpp::null ) {
      }

  };

  class selftestcase_Angle : public selftestcase {
    public:
      selftestcase_Angle(const mysqlpp::sql_int &molid, const mysqlpp::sql_int &ffid, const TopCoord& c, const TopElement &e) :
        selftestcase( 0, molid, ffid, 
            "angle", boost::lexical_cast<string>(c.f), e.i, e.j, e.k, mysqlpp::null,
            c.c0, c.c1, c.c2, mysqlpp::null, mysqlpp::null, mysqlpp::null ) {
      }

  };

  class selftestcase_Dihedral : public selftestcase {
    public:
      selftestcase_Dihedral(const mysqlpp::sql_int &molid, const mysqlpp::sql_int &ffid, const TopCoord& c, const TopElement &e) :
        selftestcase( 0, molid, ffid, 
            "dihedral", boost::lexical_cast<string>(c.f), e.i, e.j, e.k, e.l,
            c.c0, c.c1, c.c2, c.c3, c.c4, c.c5 ) {
      }

  };

  //TODO: impropers

}

#endif

