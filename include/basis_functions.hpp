#ifndef BASIS_FUNCTIONS_HPP
#define BASIS_FUNCTIONS_HPP

#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include "mesh.hpp"

template <int dim, int degree>
class basis_function {

public:

    typedef element<dim> FE;
    typedef Rd<dim> Rn;

    static const int D    = dim;        // space dimension
    static const int deg  = degree;     // polynomial degree

    int ndof;                           // number of degrees of freedom per element
    
    basis_function() {ndof = 0;}

    virtual void eval(const Rn &x, std::vector<double> &phi) const {assert(false);};

    // dphi is a matrix of size ndof x D
    virtual void eval_d(const FE &K, const Rn &x, std::vector<std::vector<double>> &dphi) const {assert(false);};

    ~basis_function() {}

};

// P1 Lagrange 1D
template <int degree>
class lagrange_1d : public basis_function<1, degree> {

    typedef typename basis_function<1, degree>::FE FE;
    typedef typename basis_function<1, degree>::Rn Rn;

public:
    
    lagrange_1d() : basis_function<1, degree>() {this->ndof = degree + 1;}

    void eval(const Rn &x, std::vector<double> &phi) const override;

    // dphi is a matrix of size ndof x D
    void eval_d(const FE &K, const Rn &x, std::vector<std::vector<double>> &dphi) const override;

    ~lagrange_1d() {}


};


#include "../src/basis_functions.tpp"


#endif 