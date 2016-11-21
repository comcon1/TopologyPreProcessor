#ifndef INCLUDE_PARAMSET_HPP_
#define INCLUDE_PARAMSET_HPP_

#include <map>

//#define TPP_UNIT_TEST

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
	typedef std::pair<std::string,std::string> Parameter; //! parameter for everything


	///
	///	\brief A dictionary of parameters. Some bicycle related to serialization.
	///
	class Parameters
	{
		std::map<std::string, std::string> _params;
	public:
		Parameters() {}
		~Parameters() {}

		template<typename T>
		bool exists(const T& key) const {  return _params.find(key) != _params.end(); }

		template<typename T>
		const std::string read(const T& key) const {
			auto ptr = _params.find(key);
			if (ptr == _params.end()) return std::string(""); //! A silent fallback. A bad idea.
			return ptr->second;
		}

		template<typename T1, typename T2>
		void add(const T1& par, const T2& val) { BOOST_CHECK( _params.insert(Parameter(par,val)).second ); }

		template<typename T>
		void remove(const T& key) {
			   BOOST_CHECK(exists(key));
			   _params.erase(_params.find(key));
			}
	}; // class Parameters


} // namespace tpp



#endif /* INCLUDE_PARAMSET_HPP_ */
