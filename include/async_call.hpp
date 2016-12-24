#ifndef H_TPP_ASYNC
#define H_TPP_ASYNC

#include <iostream>
#include <thread>
#include <future>

#include "logger.h"
#include <boost/format.hpp>

namespace tpp {

  /**
   * \brief run the function and terminate if it does not finish in the predefined time.
   */
  template<class T, typename F, class...X> T run_with_timeout(unsigned timeout, F f, X... x ) {
    auto future = std::async(std::launch::async, f, x...);
    if (future.wait_for(std::chrono::seconds(timeout)) == std::future_status::timeout) {
      TPPE << boost::format("Timeout %d failed!") % timeout; // TPP log entry will be placed here
      std::terminate();
    }
    return future.get();
  }

}

#endif
