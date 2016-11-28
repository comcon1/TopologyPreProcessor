#ifndef TPP_RUNTIME_H
#define TPP_RUNTIME_H

#include "paramset.hpp"

namespace tpp {


	/**
	 * 	\brief Some logging\ debug class, probably?
	 *
	 */
	class Runtime {
	private:
		FILE *cash;
		FILE *log;
	public:
		Runtime(const char *, const char *); // opening files log & cash
		void log_write(const char *); // string writing to log
		void log_write(std::string s) {
			log_write(s.c_str());
		}
		void cash_write(const char *, unsigned); // binary writing to cash
		~Runtime();
	};

	extern Runtime runtime; /// GLOBAL VARIABLES ARE BAD

}

#endif
