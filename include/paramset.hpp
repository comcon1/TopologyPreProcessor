#ifndef INCLUDE_PARAMSET_HPP_
#define INCLUDE_PARAMSET_HPP_

#include <map>

//
//	This whole thing should be refactored into a separate class
//


#ifdef TPP_UNIT_TEST

 #define BOOST_TEST_BUILD_INFO yes
 #define BOOST_TEST_REPORT_LEVEL detailed
 #define BOOST_TEST_REPORT_FORMAT HRF
 #define BOOST_TEST_LOG_LEVEL all

    #include <boost/test/unit_test.hpp>
    #include <boost/test/unit_test_monitor.hpp>
    #include <boost/bind.hpp>
    #include <boost/test/framework.hpp>
    #include <boost/test/unit_test_log.hpp>
    #include <boost/test/results_collector.hpp>
    #include <boost/test/results_reporter.hpp>
    #include <boost/test/test_tools.hpp>
       using boost::unit_test::test_suite;
       using boost::bind;
       using boost::ref;
       namespace ut = boost::unit_test;
 #define _BLOB_ \
                perror("UNIT_TEST is turned on so you should make all the BLOBS..\n"); \
                assert(0);

#else
 #include <cassert>
 #define _BLOB_(x) cout << "[ UNDER CONSTRUCTION - " << #x << " ]" << endl;
 // overloading of BOOST_TEST operators
 #define BOOST_CHECK(x) if (! ( x ) ) \
      printf("CHECK failed at procedure in %s at line %d\n --> %s", \
              __FILE__, __LINE__, #x);
 #define BOOST_REQUIRE(x) assert ( x );
 #define BOOST_FAIL(x) perror(#x); abort();
 #define BOOST_ERROR(x) perror(#x); abort();

 #ifdef ALLOW_WARNINGS
  #define BOOST_WARN(x) perror(#x);
 #else
  #define BOOST_WARN(x) ((void)0)
 #endif

#endif



namespace tpp
{
	// parameters for everything
	typedef std::map<std::string, std::string> t_input_params;
	typedef std::pair<std::string,std::string> t_input_param;
	// fast work with params
	template<typename T>
	bool PARAM_EXISTS(const t_input_params &pars, const T x) {
	  return (pars.find(x) != pars.end());
	}
	template<typename T>
	const std::string PARAM_READ(const t_input_params &pars, const T x) {
	  return (pars.find(x) != pars.end()) ? ((pars.find(x))->second) : std::string("");
	}
	template<typename T1, typename T2>
	void PARAM_ADD(t_input_params &pars, const T1 par, const T2 val) {
	   BOOST_CHECK( (pars.insert(t_input_param(par,val))).second );
	}
	template<typename T>
	void PARAM_DEL(t_input_params &pars, const T par) {
	   BOOST_CHECK(PARAM_EXISTS(pars,par));
	   pars.erase( pars.find(par) );
	}
}



#endif /* INCLUDE_PARAMSET_HPP_ */
