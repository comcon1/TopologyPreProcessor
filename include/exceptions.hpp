/*! \file exceptions.hpp
 */

#ifndef TPP_EXCEPTION_H
#define TPP_EXCEPTION_H

#include <string>
#include <stdexcept>
#include <map>

namespace tpp {

  /*!
   *	\brief A common exception class to be used throughout TPP project.
   *
   *	It can be copied.
   */
  class Exception : public std::exception{
    protected:
      std::string _mesg; /// basic message
      std::map<std::string, std::string> _params;
      mutable std::string _totalMessage; /// the whole message for what() will be written here
    public:

      virtual ~Exception();
      Exception(const std::string& c = "Unspecified tpp::exception");
      Exception(const char *);

      const std::map<std::string, std::string> params() const;

      /// Add a named parameter to the exception's description.
      Exception& add(const std::string param, const std::string& value);

      /// Inherited from std::runtime_error. Returns message and all parameters as one string.
      const char* what() const noexcept override;

  };

  /*!
   *	\brief A  exception class that refers to database stuff.
   *
   *  Does this class need to exists in light of exception refactor?
   *
   *	TODO: rewrite via MYSQLPP::QUERY interface.
   */
  class DbException: public Exception {
    public:
      DbException(const char *a);
      DbException(const std::string& c);
  };

}

#endif
