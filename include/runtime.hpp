#ifndef TPP_RUNTIME_H
#define TPP_RUNTIME_H

#include "paramset.hpp"

namespace tpp {

extern t_input_params cmdline; //!< what's that for?

// -- for example for command line :)

	/**
	 * 	\brief Some logging\ debug class, probably?
	 *
	 */
	class t_runtime {
	private:
		FILE *cash;
		FILE *log;
	public:
		t_runtime(const char *, const char *); // opening files log & cash
		void log_write(const char *); // string writing to log
		void log_write(std::string s) {
			log_write(s.c_str());
		}
		void cash_write(const char *, unsigned); // binary writing to cash
		~t_runtime();
	};
	extern t_runtime runtime;

}

#endif
