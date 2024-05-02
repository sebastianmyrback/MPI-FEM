#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP

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

    virtual void eval(const Rd<dim> &x, std::vector<double> &phi) const = 0;

    // dphi is a matrix of size ndof x D
    virtual void eval_d(const FE &K, const Rn &x, std::vector<std::vector<double>> &dphi) const = 0;

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

#endif 