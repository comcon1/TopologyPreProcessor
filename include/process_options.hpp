#ifndef TPP_PROCOPT_H
#define TPP_PROCOPT_H

#include <string>

namespace tpp {

  /** \brief Check if file is ready to be used as an output file.
    *
    *  @param[in,out] fname  Path to the file
    *
    *  @param[in]     ext    Extension with dot, e.g. ".txt". When specified,
    * function checks if *fname* has specified extension and change it in an
    * opposite case.
    *
    *  @param[in]     suffix Suffix to file stem (default NULL)
    */
  void processOutputWithExt(std::string &fname, const char *ext, const char *suffix = NULL);
}

#endif
