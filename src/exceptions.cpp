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
  Exception::Exception(const string& msg): _mesg(msg) {

  }

  Exception::Exception(const char* msg): _mesg(msg) {

  }

  Exception::~Exception()
  {

  }

  // The following mechanic is obsolete and should not be used.
  /*

  Exception::Exception(const char *_mesg, Parameters &_pars): mesg(_mesg), pars(_pars) {
    if (pars.exists("fatal")) {
      TPPE << format("TPP was abnormally terminated at) %s ") % boost::lexical_cast<string>(boost::posix_time::second_clock::local_time());
      TPPE << format("Position: %1% -> %2%") % pars.read("classname") % pars.read("procname");
      TPPE << mesg;
      std::terminate();
    }
  }
  */


  const char* Exception::what() const noexcept {
    _totalMessage = _mesg;
    if (_params.size() == 0){
      _totalMessage += " No parameters specified.";
    }
    else{
      _totalMessage += " Parameters:\n";
      for (const auto& param: _params){
        _totalMessage += param.first + ": " + param.second + "\n";
      }
    }
    return _totalMessage.c_str();
  }

  const std::map<std::string, std::string> Exception::params() const{
    return _params;
  }

  Exception& Exception::add(const std::string param, const std::string& value){
    _params[param] = value;
    return *this;
  }

  //
  //	DBExceptionmethods
  //
  DbException::DbException(const string& msg): Exception(msg) {

  }

  DbException::DbException(const char* msg): Exception(msg) {

  }


}
