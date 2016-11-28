#ifndef TPP_BOND_MATRIX_HEADER
#define TPP_BOND_MATRIX_HEADER

#include "core.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

namespace tpp {

/**
 * \brief Bond. Matrix Bond.
 */
class BondMatrix {
private:
  std::set<int> curset;
  boost::numeric::ublas::matrix<bool> mtx;
  void rec_check(int i);
  bool check();
public:
  /**
   * \brief Construct a bond matrix for a, mmm, something mol-related?
   *
   * \param param Surely means something.
   *
   */
  BondMatrix(const OpenBabel::OBMol& param);
  ~BondMatrix() {}

  const boost::numeric::ublas::matrix<bool> &operator *() {
    return mtx;
  }
  void print(std::ostream &os);
};
}

#endif
