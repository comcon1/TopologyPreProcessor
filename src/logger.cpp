#include "logger.hpp"

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/console.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;

namespace tpp {

  Logger logger;

  Logger& get_logger()
  {
    return logger;
  }

  void initiate_logging(const std::string& logname, bool verbose)
  {
    logging::add_file_log
      (
      keywords::file_name = logname,
      keywords::format =
      (
        expr::stream
        << "["
        << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
        << "]"
        << "[" << logging::trivial::severity
        << "]:  " << expr::smessage
        )
    );
    logging::add_console_log(std::cout,
                    boost::log::keywords::format = "%Message%",
                    boost::log::keywords::filter =
                      (!verbose && (logging::trivial::severity >= boost::log::trivial::severity_level::info) ) ||
                      (verbose &&  (logging::trivial::severity >= boost::log::trivial::severity_level::debug) )
                    );

    logging::add_common_attributes();
  }

}
