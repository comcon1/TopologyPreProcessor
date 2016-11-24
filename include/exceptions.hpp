/*! \file exceptions.hpp
 */

#ifndef TPP_EXCEPTION_H
#define TPP_EXCEPTION_H

#include "paramset.hpp"

#include <string>

namespace tpp{



	/*!
	 *	\brief A common exception class to be used throughout TPP project.
	 *
	 *	A bicycle, really. Probably better if it was inherited from std::exception.
	 */
	class Exception {
		protected:
			std::string mesg;
			Parameters pars;
	  public:
	   Exception();
	   virtual ~Exception();
	   Exception(const char *, Parameters &);
	   Exception(const char *s);
	   virtual std::string operator [] (const char *s) const;
	   virtual std::string operator [] (const std::string &s) const;

	   // rebuild all the files if you change pure-virtual function
	   virtual void fix_log() const;

	};

	/*!
	 *	\brief A  exception class that refers to database stuff.
	 *
	 *	TODO: rewrite via MYSQLPP::QUERY interface
	 */
	class DbException: public Exception {
	  public:
		DbException(const char *a, Parameters &b);
	    virtual void fix_log() const;
	};

}

#endif
