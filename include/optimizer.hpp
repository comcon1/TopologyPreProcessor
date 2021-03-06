#ifndef TPP_OPTIMIZER_H
#define TPP_OPTIMIZER_H

#include "global.hpp"

// functor suited for optimizer class
// -- functor must have an () operator
// -- with ublas::vector as argument
// -- and ublas::vector as returned value type
// -------------------------------------------


<template class Functor>
class tpp_optimizer {
  protected:
    Functor F;
    static const ublas::vector v_pred;
    static constublas::vector v_start;
  public:
    tpp_optimizer(const Functor &, const ublas::vector &, 
        const ublas::vector &, tpp_input_params &) throw (tpp_exception);
    ~tpp_optimizer();
};

#endif

    