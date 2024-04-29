#include "lagrange.hpp"

void p1_lagrange_1d::eval(const Rd<1> &x, std::vector<double> &phi) const {
    // x - quadrature point in physical element 

    phi.resize(ndof);
    phi[0] = 1.0 - x[0];
    phi[1] = x[0];
};

void p1_lagrange_1d::eval_d(const FE &K, const Rd<1> &x, std::vector<std::vector<double>> &dphi) const {
    // x - quadrature point in physical element K 
    // dphi is a matrix of size ndof x D

    dphi.resize(ndof);

    const std::vector<double> vertices = {K(0).x[0], K(1).x[0]};

    for (int i=0; i<ndof; i++) {
        int i1 = (i+1)%2;
        dphi[i][0] = 1. / (vertices[i] - vertices[i1]);
    }

};