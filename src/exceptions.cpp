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

Exception::Exception(const char *, Parameters &) {
	using std::cerr;
    if (pars.exists("fatal")) {
      cerr << "TPP was abnormally terminated at " << boost::posix_time::second_clock::local_time() << "\n";
      cerr << boost::format("Position: %1% -> %2% \n") % pars.read("classname") % pars.read("procname");
      cerr << mesg << std::endl;
      cerr << "Saving log and cache files..";
      runtime.~Runtime(); // TODO remoce explicit destructor call
      cerr << "Ok.\n";
      cerr << "---------------------------------------------------------\n";
      exit(1);
    }
}

Exception::Exception(const char *s): mesg(s){

}

std::string Exception::operator [](const char *s) const {
	return std::string(pars.read(s));
}

std::string Exception::operator [](const std::string &s) const {
	return std::string(pars.read(s));
}


void Exception::fix_log() const {
	std::ostringstream os;
	os << "TPP catched exception!\n";
	os << boost::format("***** from %1% -> %2%\n") % pars.read("classname")
					% pars.read("procname");
	os << boost::format("***** ===[ %1% ]===\n") % pars.read("error");
	if (pars.exists("line"))
		os << "***** parsing line: #" << pars.read("line") << std::endl;
	os << "***** " << mesg << std::endl;
	runtime.log_write(os.str());
}

//
//	DBExceptionmethods
//
DbException::DbException(const char *a, Parameters &b): Exception(a,b)
{

}


void DbException::fix_log() const {
	using std::endl;

	std::ostringstream os;
	os << "TPP catched exception!\n";
	os << boost::format("***** from %1% -> %2%\n")
					% pars.read("classname")
					% pars.read("procname");
	os << boost::format("***** ===[ %1% ]===\n") % pars.read("error");
	if (pars.exists("line"))
		os << "***** parsing line: #" << pars.read("line") << endl;
	os << "***** SQL error: #" << pars.read("sql_error") << endl;
	os << "***** SQL query: #" << pars.read("sql_query") << endl;
	os << "***** " << mesg << endl;
	runtime.log_write(os.str());

}
