#include "global.hpp"

namespace tpp {

t_input_params cmdline;
t_runtime runtime("log","cash");
// runtime constructor
t_runtime::t_runtime(const char *logn, const char *cashn) {
  log = fopen(logn, "a+");
  BOOST_REQUIRE(log);
  cash = fopen(cashn,"a+b");
  BOOST_REQUIRE(cash);
  cash_write("cash_start",10);
  log_write(string("\nStarting ") + PACKAGE_STRING + "\n");
  log_write(string("Build at: ") + CONFIGURE_CDATE + "\n");
  log_write(string("Log initiated at ") + lexical_cast<string>(second_clock::local_time()).c_str() + "\n");
}

// write string into log
void t_runtime::log_write(const char *s) {
  fputs(s, log);
  fflush(log);
}

// write block into cash \0 is not implemented as end-of-string
void t_runtime::cash_write(const char *s, unsigned c) {
  for (int i=0; i<c; i++)
    fputc(s[i], cash);
  fflush(cash);
}

// finishing write
t_runtime::~t_runtime() {
  log_write(string("\nLog terminated at ") + lexical_cast<string>(second_clock::local_time()) + "\nBye!");
  cash_write("cash_end",8);
  fclose(log);
  fclose(cash);
}

// exception main constructor
t_exception::t_exception(const char *s, t_input_params &p): mesg(s), pars(p) {
      if (PARAM_EXISTS(pars,"fatal")) {
        cerr << "TPP was abnormally terminated at " << second_clock::local_time() << "\n";
        cerr << format("Position: %1% -> %2% \n") % PARAM_READ(pars, "classname") % PARAM_READ(pars, "procname");
        cerr << mesg << endl;
        cerr << "Saving log and cache files..";
        runtime.~t_runtime();
        cerr << "Ok.\n";
        cerr << "---------------------------------------------------------\n";
        exit(1);
      }
}

}

