#include "db_base.hpp"

#include "global.hpp" // for magic number macro
#include "exceptions.hpp"
#include "logger.hpp"

#include <string>
#include <algorithm>

#include <boost/format.hpp>

using boost::format;

using std::string;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;
using std::to_string;
using std::ostringstream;

namespace tpp {


  /**
   * \brief Standard constructor with DB connection
   */
  DbBase::DbBase(const Settings& set) :
       settings(set), con(new mysqlpp::Connection(false)) {
    ;
  }

  DbBase::~DbBase(){
    delete con;
  }

  /**
   * \brief Standart DB connection.
   */
  bool DbBase::connectDB()  {
    TPPD << "Performing SQL connection";
    TPPD << format("Connection server: %s:%d") % settings.host % settings.port;
    TPPD << format("Authorization: %s %s") % settings.user % settings.password;
    TPPD << "Database: " + settings.dbname;
    con->connect(
        settings.dbname.c_str(),
        (settings.host + ":" + to_string(settings.port)).c_str(),
        settings.user.c_str(),
        settings.password.c_str()
        );

    // connection established
    if (!con->connected()) {
      SqlException e("SQL connection failed!");
      e.add("procname", "tpp::atom_definer::connectDB");
      e.add("error", "SQL connection error");
      e.add("sql_error", con->error() );
      throw e;
    }
    TPPD << "Connection established";

    mysqlpp::Query qu = this->con->query();
    QueryResult res;
    mysqlpp::Row row;
    ostringstream os,os2;

    TPPD << "Checking DB version." << endl;
    os2 << "SELECT `id`,`keyword`,`value` FROM `properties` WHERE keyword='magic_number'";
    #ifdef SQLDEBUG
    TPPD << os2.str();
    #endif // SQLDEBUG
    qu << os2.str();
    res = qu.store();
    if (!res) {
      SqlException e("SQL query failed: may be you have incorrect DataBase!");
      e.add("procname", "tpp::DbInfo::connectDB");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error() );
      e.add("query", os2.str());
      throw e;
    }
    if (res.num_rows() != 1) {
      Exception e("Your DataBase is empty!");
      e.add("procname", "tpp::DbInfo::connectDB");
      e.add("error", "Error in DB check.");
      throw e;
    }

    string cur_mn(res.at(0)["value"]);
    string required_mn(DB_MAGICNUMBER);

    os << format("Required|current magic number: %s | %s") % required_mn % cur_mn;

    TPPD << os.str();

    if (cur_mn.substr(0,4) != required_mn.substr(0,4)) {
      Exception e("Your DataBase has incompatible version!");
      e.add("procname", "tpp::DbInfo::connectDB");
      e.add("error", "Error in DB check.");
      throw e;
    }

    if (cur_mn.c_str()[4] > required_mn.c_str()[4]) {
      TPPE << "\n"
"** Your DataBase version is slightly higher. Program should work but\n"
"** this behavior may be unpredictable. It is better to update the code.\n";
    }

    return true;
  } // end connectDB

  /**
    * \brief DB INFO implementation
    */
  DbInfo::DbInfo(const Settings& sets, const std::string& ffn): DbBase(sets), ffName(ffn) {
    ;
  }

  /**
   * \brief DB queries for DB-INFO class
   */
  bool DbInfo::connectDB() {
    DbBase::connectDB();
    this->loadFFData();
    return true;
  }

  /** \brief Load FF data from DB and set up object variables
    */
  void DbInfo::loadFFData() {
    mysqlpp::Query qu = this->con->query();
    QueryResult res;
    mysqlpp::Row row;
    ostringstream os;
    // get ff info
    os << format("SELECT `id`,`include`,`desc`, `nbfunc`, `combrule`, "
    "`genpairs`, `fudgelj`, `fudgeqq` FROM `forcefield` WHERE name='%1$s'") % this->ffName.c_str();
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    qu << os.str();
    res = qu.store();
    if (!res) {
      SqlException e("SQL query failed!");
      e.add("procname", "tpp::DbInfo::loadFFData");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error() );
      e.add("query", qu.str());
      throw e;
    }
    if (!res.num_rows()) {
      Exception e("Force field  not found!");
      e.add("procname", "tpp::DbInfo::loadFFData");
      e.add("error", "Error in parameters");
      e.add("field", this->ffName);
      throw e;
    }
    this->ffID = res.at(0)["id"];
    this->ffInclude = (res.at(0)["include"]).c_str();
    this->ffDesc = (res.at(0)["desc"]).c_str();

    os.str("");
    os << format("%1d %1d %s %5.4f %5.4f") % (int) (res.at(0)["nbfunc"]) % (int) (res.at(0)["combrule"])
      % ( ( (bool) (res.at(0)["genpairs"]) ) ? "yes" : "no" )
      %  (double) (res.at(0)["fudgelj"]) % (double) (res.at(0)["fudgeqq"]);
    this->ffDefaults = os.str();

    qu.reset();
    os.str("");
    os << format("select `value` from `properties` WHERE `keyword`='%1$srev'") % this->ffName;
    qu << os.str();
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    res = qu.store();
    if ( (!res) || res.size() == 0 || (!res.at(0)) ) {
      SqlException e("Error in SQL or force field revision is unknown!");
      e.add("procname", "tpp::atom_definer::connectDB");
      e.add("error", (!res) ? "SQL query error" : "Revision not found");
      e.add("sql_error", qu.error() );
      e.add("query", os.str());
      throw e;
    }
    this->ffRev = res.at(0)["value"].c_str();
  } // end loadFFData


  string DbInfo::getStatistics() {

    mysqlpp::Query qu = this->con->query();
    QueryResult res;
    mysqlpp::Row row;
    ostringstream os, retos;
    os << format("\
SELECT\n\
 (SELECT COUNT(*) FROM atoms WHERE ffield = %1$d) as count_atoms,\n\
 (SELECT COUNT(*) FROM bonds WHERE ffield = %1$d) as count_bonds,\n\
 (SELECT COUNT(*) FROM angles WHERE ffield = %1$d) as count_angles,\n\
 (SELECT COUNT(*) FROM dihedrals WHERE ffield = %1$d) as count_dihedrals,\n\
 (SELECT COUNT(*) as CF FROM nonbonded WHERE ffield = %1$d) as count_nonbond;") % this->ffID;
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    qu << os.str();
    res = qu.store();
    if ( (!res) || (!res.at(0)) ) {
      SqlException e( "SQL query failed!");
      e.add("procname", "tpp::atom_definer::connectDB");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error() );
      e.add("query", os.str());
      throw e;
    }
    retos << format("\
Forcefield %1$s was found in database.\n\
Description: %7$s.\n\
Total statistics:\n\
%2$5d atoms,     %3$5d bonds, %4$5d angles,\n\
%5$5d dihedrals, %6$5d nonbonded parameters.\n")
      % this->ffName.c_str() % res.at(0)["count_atoms"] % res.at(0)["count_bonds"]
      % res.at(0)["count_angles"] % res.at(0)["count_dihedrals"] % res.at(0)["count_nonbond"]
      % this->ffDesc;
    //TODO: this query does not show last INSERT events. Should fix.
    qu.reset();
    os.str("");
    os << format("show table status from `%1$s`") % settings.dbname;
    #ifdef SQLDEBUG
    TPPD << os.str();
    #endif // SQLDEBUG
    qu << os.str();
    res = qu.store();
    if ( (!res) || (!res.at(0)) ) {
      SqlException e("SQL status query failed!");
      e.add("procname", "tpp::atom_definer::connectDB");
      e.add("error", "SQL query error");
      e.add("sql_error", qu.error() );
      e.add("query", os.str());
      throw e;
    }

    vector<string> ss;
    for(mysqlpp::Row::size_type co = 0; co < res.num_rows(); ++co) {
      row = res.at(co);
      if (row["Update_time"] != "NULL")
        ss.push_back(string(row["Update_time"]));
    }
    std::sort(ss.begin(), ss.end());
    retos << "Database last update: " << ss.back() << endl;

    // force field revision
    retos << format("Force field %1$s DB revision: %2$s.") % this->ffName % this->ffRev;
    return retos.str();
  } // end get Statistics

}
