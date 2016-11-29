/** \file db_scanner.hpp
 *
 *	\brief One more mysterious header that has something to do with parsing...
 *
 */

#ifndef TPP_DBSCANNER_H
#define TPP_DBSCANNER_H

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

  /// for brevity
  typedef mysqlpp::StoreQueryResult QueryResult;

  /**
   *
   *  \brief Class that encapsulates basic working with databases.
   *
   */
  class DbBase {

    public:
      /// Database connection settings.
      struct Settings
      {
        std::string host;
        std::string user;
        std::string password;
        unsigned port;
        std::string dbname;
      };

      DbBase(const Settings& settings);
      virtual ~DbBase();

    protected:
      Settings settings;
      mysqlpp::Connection *con;
      virtual bool connect_db();
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
      DbInfo(const Settings& set, const std::string& ffn);
      int get_ffid() { return ffid; }
      std::string get_ffinclude() { return ffinclude; }
      std::string get_ffrev() { return ffrev; }

      //! Function that gives string statistics about current force field.
      std::string get_statistics();

  };

}

#endif

