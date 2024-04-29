#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP

#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include "mesh.hpp"

// P1 Lagrange 1D
class p1_lagrange_1d {
    typedef element<1> FE;

public:
    static const int D    = 1;
    static const int ndof = 2;
    static const int deg  = 1;


    p1_lagrange_1d() {}

    void eval(const Rd<1> &x, std::vector<double> &phi) const;

    // dphi is a matrix of size ndof x D
    void eval_d(const FE &K, const Rd<1> &x, std::vector<std::vector<double>> &dphi) const;

    ~p1_lagrange_1d() {}


};

#endif 