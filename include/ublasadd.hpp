// http://osdir.com/ml/lib.boost.ublas/2005-12/msg00016.html
// from boost::ublas mailing list archive (2005 year)
// *********************************************************
#ifndef UBLAS_ADD_H
#define UBLAS_ADD_H

#include "boost/numeric/ublas/vector.hpp" 
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/triangular.hpp"
#include "boost/numeric/ublas/operation.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/bindings/lapack/gesvd.hpp"
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"

template<class M> 
 int svd_rank(M const& m) {
   namespace ublas = boost::numeric::ublas;
   namespace lapack = boost::numeric::bindings::lapack;
   if (m.size2() == 1)
     for (int i=0; i<m.size1(); i++)
       if (fabs(m(i,0)) > 1e-10) return 1;
   if (m.size1() == 1)
     for (int i=0; i<m.size2(); i++)
       if (fabs(m(0,i)) > 1e-10) return 1;
   ublas::matrix<double,ublas::column_major> _m(m); // create working copy
   ublas::vector<double> _v(m.size1() < m.size2() ? m.size1() : m.size2());
   lapack::gesvd(_m, _v);
   int rank = 0;
   for (int i=0; i<_v.size(); ++i)
     if (fabs(_v(i)) > 1e-10) rank++;
   return rank;
 }

template<class M1, class M2>
void svd_inv(M1 const& m, M2& inv) {
   BOOST_REQUIRE(m.size1() == m.size2());
   namespace ublas = boost::numeric::ublas;
   namespace lapack = boost::numeric::bindings::lapack;  
   ublas::matrix<double, ublas::column_major> _m(m), _v(m.size1(), m.size2()), _u(m.size1(), m.size2());
   ublas::vector<double> s00(m.size1());
   lapack::gesvd(_m, s00, _u, _v);
   _v = trans(_v);
   _u = trans(_u);
   ublas::matrix<double, ublas::column_major> S00(m.size1(), m.size2());
   S00 *= 0.0;
   for (int i=0; i<s00.size(); ++i) S00(i,i) = 1/s00(i);
   inv = prod(_v,S00);
   inv = prod(inv, _u);
   return;
}
// determinant
template<class M>
 double lu_det(M const& m) {
   using namespace boost::numeric::ublas;
   // create a working copy of the input
   matrix<double> mLu(m);
   permutation_matrix<std::size_t> pivots(m.size1());
   lu_factorize(mLu, pivots);

   double det = 1.0;

   for (std::size_t i=0; i < pivots.size(); ++i) {
     if (pivots(i) != i)
       det *= -1.0;
     det *= mLu(i,i);
   }
   return det;
 }

template<class M1, class M2>
void lu_inv(M1 const& m, M2& inv) {
  using namespace boost::numeric::ublas;
  BOOST_REQUIRE(m.size1() == m.size2() );
  BOOST_REQUIRE( (inv.size1() == m.size1()) && (inv.size1() == m.size2()) );

  // create a working copy of the input
  matrix<double> mLu(m);

  // perform LU-factorization
  lu_factorize(mLu);

  // create identity matrix of "inverse"
  inv.assign(identity_matrix<double>(m.size1()));

  // backsubstitute to get the inverse
  lu_substitute<matrix<double> const, M2 >(mLu, inv);
}

#endif
