#include <iostream>
#include "../include/strutil.hpp"
#include "boost/filesystem.hpp"

using std::string;
using std::cout;
using std::endl;
namespace bfs = boost::filesystem;

int main() {
   string file = "../config.11";
   auto vec = strutil::split(file,".");
   string v = vec.back();
   cout << "LAST: " << v << endl;
   bfs::path p (file);
   cout << "BOOST FS: " << p.extension().string() << endl;
   cout << "BOOST BN: " << p.stem().string() << endl;
   cout << "CURRENT: " << bfs::current_path() << endl;
   cout << bfs::weakly_canonical(p) << endl;
   p = bfs::path(p.string() + ".rtp");
   cout << bfs::weakly_canonical(p) << endl;
   return 0;
}
