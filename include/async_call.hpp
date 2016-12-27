#ifndef H_TPP_ASYNC
#define H_TPP_ASYNC

#include <iostream>
#include <thread>
#include <future>
#include <ctime>

#include "logger.hpp"
#include <boost/format.hpp>

namespace tpp {

  /**
   * \brief run the function and terminate if it does not finish in the predefined time.
   */
  template<class T, typename F, class...X> T run_with_timeout(unsigned timeout, F f, X... x ) {
    TPPD << " << Starting asynchronous call >>";
    std::clock_t beg = std::clock();
    auto future = std::async(std::launch::async, f, x...);
    if (future.wait_for(std::chrono::seconds(timeout)) == std::future_status::timeout) {
      TPPE << boost::format(" << Timeout %d failed! >>") % timeout; // TPP log entry will be placed here
      std::terminate();
    }
    std::clock_t en = std::clock();
    TPPD << boost::format(" << Asynchronous call finished in %.3f sec. >>") % ((double)(en - beg) / CLOCKS_PER_SEC);
    return future.get();
  }

  // TODO: envelop for asynchronous run of class methods
  // now we can use C++11 lambda instead like this:
  //
  // std::cout << run_with_timeout<int>(3, 
  //          [&]() { std::cout << "**** In da lambda!!" << std::endl; return f(100); }
  //           ) << std::endl;


}

#endif
