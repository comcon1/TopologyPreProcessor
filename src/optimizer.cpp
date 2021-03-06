#include "optimizer.hpp"

tpp_optimizer::tpp_optimizer(Functor f, ublas::vector v0, 
    ublas::vector v1) throw (tpp_exception):  F(f), v_start(v0), v_pred(v1) {
  // example call
  if (F(v0).size() != v1.size()) 
    throw tpp_exception("tpp_optimizer: dimension of functor returned vector doesn't match start vector.");
  
}

       