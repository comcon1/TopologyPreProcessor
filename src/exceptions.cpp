#include "exceptions.hpp"
#include "runtime.hpp"

#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace tpp;



//
// Exception methods
//
Exception::Exception(): mesg("Undefined exception.") {

}

Exception::~Exception()
{

}

Exception::Exception(const char *, t_input_params &) {
	using std::cerr;
    if (PARAM_EXISTS(pars,"fatal")) {
      cerr << "TPP was abnormally terminated at " << boost::posix_time::second_clock::local_time() << "\n";
      cerr << boost::format("Position: %1% -> %2% \n") % PARAM_READ(pars, "classname") % PARAM_READ(pars, "procname");
      cerr << mesg << std::endl;
      cerr << "Saving log and cache files..";
      runtime.~t_runtime();
      cerr << "Ok.\n";
      cerr << "---------------------------------------------------------\n";
      exit(1);
    }
}

Exception::Exception(const char *s): mesg(s){

}

std::string Exception::operator [](const char *s) const {
	return std::string(PARAM_READ(pars, s));
}

std::string Exception::operator [](const std::string &s) const {
	return std::string(PARAM_READ(pars, s));
}


void Exception::fix_log() const {
	std::ostringstream os;
	os << "TPP catched exception!\n";
	os << boost::format("***** from %1% -> %2%\n") % PARAM_READ(pars, "classname")
					% PARAM_READ(pars, "procname");
	os << boost::format("***** ===[ %1% ]===\n") % PARAM_READ(pars, "error");
	if (PARAM_EXISTS(pars, "line"))
		os << "***** parsing line: #" << PARAM_READ(pars, "line") << std::endl;
	os << "***** " << mesg << std::endl;
	runtime.log_write(os.str());
}

//
//	DBExceptionmethods
//
DbException::DbException(const char *a, t_input_params &b): Exception(a,b)
{

}


void DbException::fix_log() const {
	using std::endl;

	std::ostringstream os;
	os << "TPP catched exception!\n";
	os << boost::format("***** from %1% -> %2%\n")
					% PARAM_READ(pars, "classname")
					% PARAM_READ(pars, "procname");
	os << boost::format("***** ===[ %1% ]===\n") % PARAM_READ(pars, "error");
	if (PARAM_EXISTS(pars, "line"))
		os << "***** parsing line: #" << PARAM_READ(pars, "line") << endl;
	os << "***** SQL error: #" << PARAM_READ(pars, "sql_error") << endl;
	os << "***** SQL query: #" << PARAM_READ(pars, "sql_query") << endl;
	os << "***** " << mesg << endl;
	runtime.log_write(os.str());

}
