#include "core.hpp"
#include "logger.hpp"
#include "process_options.hpp"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <string>
#include "strutil.hpp"

namespace tpp {

  namespace bfs = boost::filesystem;
  using std::string;
  using boost::format;

  void processOutputWithExt(string &fname, const char *ext, const char* suffix) {
      bfs::path ro_path(fname);
      bfs::path ro_parent(ro_path.parent_path());
      ro_path = bfs::weakly_canonical(ro_path);
      string ro_ext = strutil::toLower( ro_path.extension().string() );
      ro_parent = ro_path.parent_path();
      if ( (ro_ext != ext) && ext) ro_ext = ext;
      if (suffix) ro_ext = string(suffix) + ro_ext;
      ro_parent /= ro_path.stem().string()+ro_ext;
      string upro = strutil::toUpper(ro_ext.substr(1));
      TPPI << format("Output %s: %s") % strutil::toUpper(ro_ext)
              % bfs::absolute(ro_parent).string();
      if ( bfs::exists(ro_path) ) {
        TPPI << "Output file exists and will be overwritten." ;
      }
  }

  void processInput(string &fname) {
  }

}
