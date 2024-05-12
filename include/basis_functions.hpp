#pragma once


#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include "mesh.hpp"

template <int dim, int degree>
class BasisFunction {

public:

    typedef Cell<dim> FE;
    typedef Point<dim> Rn;

    static const int D    = dim;        // space dimension
    static const int deg  = degree;     // polynomial degree

    int ndof;                           // number of degrees of freedom per element
    
    BasisFunction() {ndof = 0;}

    virtual void eval(const Rn &x, std::vector<double> &phi) const {assert(false);};

    // dphi is a matrix of size ndof x D
    virtual void eval_d(const FE &K, const Rn &x, std::vector<std::vector<double>> &dphi) const {assert(false);};

    ~BasisFunction() {}

};

// P1 Lagrange 1D
template <int degree>
class P1Lagrange1D : public BasisFunction<1, degree> {

    typedef typename BasisFunction<1, degree>::FE FE;
    typedef typename BasisFunction<1, degree>::Rn Rn;

public:
    
    P1Lagrange1D() : BasisFunction<1, degree>() {this->ndof = degree + 1;}

    void eval(const Rn &x, std::vector<double> &phi) const override;

    // dphi is a matrix of size ndof x D
    void eval_d(const FE &K, const Rn &x, std::vector<std::vector<double>> &dphi) const override;

    ~P1Lagrange1D() {}


};


#include "basis_functions.tpp"

