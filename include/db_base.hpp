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

namespace tpp {

  //! Exception to be used when something goes wrong on sql level.
  class SqlException : public Exception {
    public:
    SqlException(const char *s): Exception(s) { }
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

      /** \brief Initialize class with settings of DB connection
        *
        */
      DbBase(const Settings& settings);

      /** \brief Safety disconnect DB
        */
      virtual ~DbBase();

      /** \brief Connect to database.
        *  Method is called from the constructor.
        */
      virtual bool connectDB();

    protected:
      Settings settings;         //!< Structure with DB connection settings
      mysqlpp::Connection *con;  //!< Pointer to Connection object

  };

  /**
   * \brief Class that accepts information about DB and Force Field (FF).
   */
  class DbInfo: public DbBase {

    protected:

      int ffID;                   //!< ID of the force field in DB.
      std::string ffName;         //!< Name of the force field
      std::string ffDesc;         //!< USE SOMEWHERE?? Description of the force field (add to ITP header)
      std::string ffInclude;      //!< DEPRECATED
      std::string ffDefaults;     //!< defaults string for forcefield.itp

      /** \brief Revision of the force field
        *
        * Very important force field characteristic. One should be sure that force field
        * is the same in data base and in the repository.
        */
      std::string ffRev;

      /** \brief Load FF data from DB and set up object variables
        */
      void loadFFData();

    public:
      DbInfo(const Settings& set, const std::string& ffn);

      /** \brief Overriden method of DB connection.
        *
        * Loads FF information via loadFFData call.
        */
      bool connectDB() override;
      
      int getFFID() { return ffID; }                    //!< Public alias to ffID
      std::string getFFInclude() { return ffInclude; }  //!< Public alias to ffInclude [DEPRECATED]
      std::string getFFRev() { return ffRev; }          //!< Public alias to ffRev
      std::string getFFDefaults() { return ffDefaults; }//!< String of FF defaults [for EXPANDED mode]

      /** \brief Function that gives string statistics about current force field.
        *
        * This function performs call curious query that may be sensitive
        * to SQL version
        */
      std::string getStatistics();

  };

}

#endif

