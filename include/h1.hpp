#ifndef LAGRANGE_HPP
#define LAGRANGE_HPP

#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include "mesh.hpp"

class h1_1d {
    typedef element<1> FE;

public:

    static const int ndof = 2;
    static const int deg = 1;

    h1_1d() {}

    void eval(const double &x, std::vector<double> &phi) const;

    void eval_d(const double &x, std::vector<double> &dphi) const;

    ~h1_1d() {}


};

#endif 