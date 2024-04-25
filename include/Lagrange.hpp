#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP

#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include "mesh.hpp"

class P1Lagrange1D {
    typedef element<1> FE;

public:

    static const int ndof = 2;
    static const int deg = 1;

    P1Lagrange1D() {}

    void eval(const double &x, std::vector<double> &phi) const;

    void eval_d(const double &x, std::vector<double> &dphi) const;

    ~P1Lagrange1D() {}


};

#endif 