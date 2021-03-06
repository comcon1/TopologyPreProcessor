#ifndef BABELRD_TPP_H
#define BABELRD_TPP_H

#include "global.hpp"

#ifndef BABEL_DIR
 #error Babel share directory not defined!
#endif

namespace tpp {
  class BabelReader {
    private:
      const curdir = BABEL_DIR;
      map<string,string> 
    public:
      BabelReader();
      ~BabelReader();
  };
}
#endif
     