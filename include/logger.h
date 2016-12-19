#pragma once

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>

#define TPPI BOOST_LOG_SEV(tpp::get_logger(), boost::log::trivial::severity_level::info)
#define TPPD BOOST_LOG_SEV(tpp::get_logger(), boost::log::trivial::severity_level::debug)
#define TPPE BOOST_LOG_SEV(tpp::get_logger(), boost::log::trivial::severity_level::error)

namespace tpp {

 typedef boost::log::sources::severity_logger<boost::log::trivial::severity_level> Logger;
//
//	This functions configures boost::log to create two output
// sources (sinks): first is written into provided file, second
// is wriiten to stdout. All messages are written into file
// (with timestamps and severity level) and messages with "info"
// or above level is written to stdout (withput timestamps and 
// severity levels)
//
 void initiate_logging(const std::string& log_file_name);

 Logger& get_logger();
}
