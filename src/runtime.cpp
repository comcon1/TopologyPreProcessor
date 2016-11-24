#include "global.hpp"
#include "runtime.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using std::string;

namespace tpp {

// Global varibales. THEY ARE BAD
Parameters cmdline;
Runtime runtime("log", "cash");
//


// runtime constructor
Runtime::Runtime(const char *logn, const char *cashn) {
	log = fopen(logn, "a+");
	BOOST_REQUIRE(log);
	cash = fopen(cashn, "a+b");
	BOOST_REQUIRE(cash);
	cash_write("cash_start", 10);
	log_write(string("\nStarting ") + PACKAGE_STRING + "\n");
	log_write(string("Build at: ") + CONFIGURE_CDATE + "\n");
	log_write(
			string("Log initiated at ")
					+ boost::lexical_cast<string>(boost::posix_time::second_clock::local_time()).c_str()
					+ "\n");
}

// write string into log
void Runtime::log_write(const char *s) {
	fputs(s, log);
	fflush(log);
}

// write block into cash \0 is not implemented as end-of-string
void Runtime::cash_write(const char *s, unsigned c) {
	for (int i = 0; i < c; i++)
		fputc(s[i], cash);
	fflush(cash);
}

// finishing write
Runtime::~Runtime() {
	log_write(
			string("\nLog terminated at ")
					+ boost::lexical_cast<string>(boost::posix_time::second_clock::local_time())
					+ "\nBye!");
	cash_write("cash_end", 8);
	fclose(log);
	fclose(cash);
}

}

