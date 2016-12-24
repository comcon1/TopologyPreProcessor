#include "exceptions.hpp"
#include "logger.hpp"

#include <exception>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>

using std::endl;
using boost::format;
using std::string;

namespace tpp {

  //
  Exception::Exception(): mesg("Undefined exception.") {

  }

  Exception::~Exception()
  {

  }

  Exception::Exception(const char *_mesg, Parameters &_pars): mesg(_mesg), pars(_pars) {
    if (pars.exists("fatal")) {
      TPPE << format("TPP was abnormally terminated at) %s ") % boost::lexical_cast<string>(boost::posix_time::second_clock::local_time());
      TPPE << format("Position: %1% -> %2%") % pars.read("classname") % pars.read("procname");
      TPPE << mesg;
      std::terminate();
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
          os << format("***** from %1% -> %2%\n") % pars.read("classname")
                                          % pars.read("procname");
          os << format("***** ===[ %1% ]===\n") % pars.read("error");
          if (pars.exists("line"))
                  os << "***** parsing line: #" << pars.read("line") << std::endl;
          os << "***** " << mesg << std::endl;
          TPPE << os.str();
  }

  //
  //	DBExceptionmethods
  //
  DbException::DbException(const char *a, Parameters &b): Exception(a,b)
  {

  }


  void DbException::fix_log() const {
          std::ostringstream os;
          os << "TPP catched exception!\n";
          os << format("***** from %1% -> %2%\n")
                                          % pars.read("classname")
                                          % pars.read("procname");
          os << format("***** ===[ %1% ]===\n") % pars.read("error");
          if (pars.exists("line"))
                  os << "***** parsing line: #" << pars.read("line") << endl;
          os << "***** SQL error: #" << pars.read("sql_error") << endl;
          os << "***** SQL query: #" << pars.read("sql_query") << endl;
          os << "***** " << mesg << endl;
          TPPE << os.str();
  }

}
