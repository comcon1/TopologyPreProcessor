#include "db_scanner.hpp"
#include "exceptions.hpp"
#include "runtime.hpp"

#include <algorithm>

#include <boost/format.hpp>

using boost::format;

using std::string;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;
using std::ostringstream;


namespace tpp {

  /*
   * SQL Exception message format.
   */
  void t_sql_exception::fix_log() const {
    std::ostringstream os; 
    os << "TPP catched exception!\n";
    os << format("***** from %1% -> %2%\n") % PARAM_READ(pars, "classname") % PARAM_READ(pars, "procname");
    os << "***** " << mesg << endl;
    os << "***** MYSQL: " << PARAM_READ(pars, "sql_error") << endl;
    runtime.log_write(os.str());
  }

  /*
   * Standard constructor with DB connection
   */
  db_base::db_base(t_input_params p_) throw (Exception):
       par(p_), con(new mysqlpp::Connection(false)) {
         ;
  }

  /*
   * Standart DB connection.
   */
  bool db_base::connect_db() throw (Exception) {
    con->connect(
        PARAM_READ(par,"dbname").c_str(),
        (PARAM_READ(par,"host")+string(":")+PARAM_READ(par,"port")).c_str(),
        PARAM_READ(par,"user").c_str(),
        PARAM_READ(par,"password").c_str()
        );

    // connection established
    if (!con->connected()) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::atom_definer::connect_db");
      PARAM_ADD(params, "error", "SQL connection error");
      PARAM_ADD(params, "sql_error", con->error() );
      throw t_sql_exception("SQL connection failed!", params);
    }

    mysqlpp::Query qu = this->con->query();
    MYSQLPP_RESULT res;
    mysqlpp::Row row;
    ostringstream os; 

    os << "Checking DB version." << endl;

    qu << "SELECT `id`,`keyword`,`value` FROM `properties` WHERE keyword='magic_number'";
    res = qu.store();
    if (!res) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::db_info::connect_db");
      PARAM_ADD(params, "error", "SQL query error");
      PARAM_ADD(params, "sql_error", qu.error() );
      throw t_sql_exception("SQL query failed: may be you have incorrect DataBase!", params);
    }
    if (res.num_rows() != 1) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::db_info::connect_db");
      PARAM_ADD(params, "error", "Error in DB check.");
      throw Exception("Your DataBase is empty!", params);
    }

    string cur_mn(res.at(0)["value"]);
    string required_mn(DB_MAGICNUMBER);
    BOOST_CHECK(cur_mn == required_mn);

    os << "Required magic number: " << required_mn << endl;
    os << "Current  magic number: " << cur_mn << endl;

    runtime.log_write(os.str()); 

    if (cur_mn.c_str()[4] > required_mn.c_str()[4]) {
      cerr << "\n"
"Your DataBase version is slightly higher. Program should work but\n"
"this behavior may be unpredictable. It is better to update the code." << endl;
    }

    return true;
  }

  /*
   * Initializing DB-INFO class.
   */
  db_info::db_info(t_input_params p) throw (Exception): db_base(p) {
    this->connect_db();
  }

  /*
   * DB queries for DB-INFO class
   */
  bool db_info::connect_db() throw (Exception) {
    db_base::connect_db();
    this->getFFdata();
  }

  /*
   * Protected method gathering force field information.
   */
  void db_info::getFFdata() {
    mysqlpp::Query qu = this->con->query();
    MYSQLPP_RESULT res;
    mysqlpp::Row row;
    // get ffid
    this->ffname = PARAM_READ(par,"ffname");
    qu << format("SELECT `id`,`include`,`desc` FROM `forcefield` WHERE name='%1$s'") % this->ffname.c_str();
    res = qu.store();
    if (!res) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::db_info::getFFdata");
      PARAM_ADD(params, "error", "SQL query error");
      PARAM_ADD(params, "sql_error", qu.error() );
      throw t_sql_exception("SQL query failed!", params);
    }
    if (!res.num_rows()) {
      t_input_params params;
      PARAM_ADD(params, "procname", "tpp::db_info::getFFdata");
      PARAM_ADD(params, "error", "Error in parameters");
      throw Exception((string("Force field '")+this->ffname+string("' not found!")).c_str(), params);
    }
    this->ffid = res.at(0)["id"];
    this->ffinclude = (res.at(0)["include"]).c_str();
    this->ffdesc = (res.at(0)["desc"]).c_str();

    qu.reset();
    qu << format("select `value` from `properties` WHERE `keyword`='%1$srev'") % this->ffname;
    res = qu.store();
    if ( (!res) || (!res.at(0)) ) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::atom_definer::connect_db");
        PARAM_ADD(params, "error", "SQL query error");
        PARAM_ADD(params, "sql_error", qu.error() );
        throw t_sql_exception("SQL force field revision is unknown!", params);
    }
    this->ffrev = res.at(0)["value"].c_str();
  }

  /*
   * Public function giving string statistics about current force field.
   */
  string db_info::get_statistics() {

    mysqlpp::Query qu = this->con->query();
    MYSQLPP_RESULT res;
    mysqlpp::Row row;

    qu << format("\
SELECT\n\
 (SELECT COUNT(*) FROM atoms WHERE ffield = %1$d) as count_atoms,\n\
 (SELECT COUNT(*) FROM bonds WHERE ffield = %1$d) as count_bonds,\n\
 (SELECT COUNT(*) FROM angles WHERE ffield = %1$d) as count_angles,\n\
 (SELECT COUNT(*) FROM dihedrals WHERE ffield = %1$d) as count_dihedrals,\n\
 (SELECT COUNT(*) as CF FROM nonbonded WHERE ffield = %1$d) as count_nonbond;") % this->ffid;
    res = qu.store();
    if ( (!res) || (!res.at(0)) ) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::atom_definer::connect_db");
        PARAM_ADD(params, "error", "SQL query error");
        PARAM_ADD(params, "sql_error", qu.error() );
        throw t_sql_exception("SQL query failed!", params);
    }
    std::ostringstream os;
    os << format("\
Forcefield %1$s was found in database.\n\
Description: %7$s.\n\
Total statistics:\n\
%2$5d atoms,     %3$5d bonds, %4$5d angles,\n\
%5$5d dihedrals, %6$5d nonbonded parameters.\n")
      % this->ffname.c_str() % res.at(0)["count_atoms"] % res.at(0)["count_bonds"]
      % res.at(0)["count_angles"] % res.at(0)["count_dihedrals"] % res.at(0)["count_nonbond"]
      % this->ffdesc;

    //TODO: this query does not show last INSERT events. Should fix.
    qu.reset();
    qu << format("show table status from `%1$s`") % PARAM_READ(par,"dbname");
    res = qu.store();
    if ( (!res) || (!res.at(0)) ) {
        t_input_params params;
        PARAM_ADD(params, "procname", "tpp::atom_definer::connect_db");
        PARAM_ADD(params, "error", "SQL query error");
        PARAM_ADD(params, "sql_error", qu.error() );
        throw t_sql_exception("SQL status query failed!", params);
    }

    vector<string> ss;
    for(mysqlpp::Row::size_type co = 0; co < res.num_rows(); ++co) {
      row = res.at(co);
      if (row["Update_time"] != "NULL")
        ss.push_back(string(row["Update_time"]));
    }
    std::sort(ss.begin(), ss.end());
    os << "Database last update: " << ss.back() << endl;

    // force field revision
    os << format("Force field %1$s DB revision: %2$s.") % this->ffname % this->ffrev  << endl;
    return os.str();
  }

}
