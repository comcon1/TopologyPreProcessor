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
  class SqlException : public Exception {
    public:
    virtual void fix_log() const; 
    SqlException(const char *s, Parameters &p): Exception(s, p) { ; }
  };

  /**
   *
   *  \brief Class that encapsulates basic working with databases.
   *
   */
  class DbBase {
    protected:
      mysqlpp::Connection *con;
      Parameters par;
      virtual bool connect_db();
    public:
      // need parameters 'host','user','dbname','password','port','ffname'
      DbBase(Parameters p);
      virtual ~DbBase() { delete con; }
  };
 
  /**
   * \brief Class that accepts information about DB and Force Field (FF).
   */
  class DbInfo: public DbBase {

    protected:
      bool connect_db() override;
      int ffid;
      std::string ffname;
      std::string ffdesc;
      std::string ffinclude;
      std::string ffrev;

      //! Protected method for gathering force field information.
      void getFFdata();
      void getDBdata();

    public:
      DbInfo(Parameters);
      int get_ffid() { return ffid; }
      std::string get_ffinclude() { return ffinclude; }
      std::string get_ffrev() { return ffrev; }

      //! Function that gives string statistics about current force field.
      std::string get_statistics();

  };

}

#endif

