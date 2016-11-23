/** \file db_scanner.hpp
 *
 *	\brief One more mysterious header that has something to do with parsing...
 *
 */

#ifndef TPP_DBSCANNER_H
#define TPP_DBSCANNER_H

#include "global.hpp"
#include "core.hpp"
#include "exceptions.hpp"

#include <mysql++/mysql++.h>
#include <set>
#include <map>


#define TPP_SMART_COEFF 10
#define TPP_ZNUC_COEFF  10000
#define TPP_BOND_COEFF  2
#define TPP_ANGLE_COEFF 4
#define TPP_DIHED_COEFF 8

namespace tpp {

  //! Exception to be used when something goes wrong on sql level.
  class t_sql_exception : public Exception {
    public:
    virtual void fix_log() const; 
    t_sql_exception(const char *s, Parameters &p): Exception(s, p) { ; }
  };

  /**
   *
   *  \brief Class that encapsulates working with databases.
   *
   */
  class db_base {
    protected:
      mysqlpp::Connection *con;
      Parameters par;
      virtual bool connect_db();
    public:
      // need parameters 'host','user','dbname','password','port','ffname'
      db_base(Parameters p);
      virtual ~db_base() { delete con; }
  };
 
  /**
   * \brief Class that accepts information about DB and FF.
   */
  class db_info: public db_base {

    protected:
      virtual bool connect_db();
      int ffid;
      std::string ffname;
      std::string ffdesc;
      std::string ffinclude;
      std::string ffrev;

      void getFFdata();
      void getDBdata();

    public:
      db_info(Parameters);
      int get_ffid() { return ffid; }
      std::string get_ffinclude() { return ffinclude; }
      std::string get_ffrev() { return ffrev; }
      std::string get_statistics();

  };

}

#endif

