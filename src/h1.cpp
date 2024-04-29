#include "h1.hpp"

void h1_1d::eval(const double &x, std::vector<double> &phi) const {
    phi.resize(ndof);
    phi[0] = 1.0 - x;
    phi[1] = x;
};

void h1_1d::eval_d(const double &x, std::vector<double> &dphi) const {
    
    dphi.resize(ndof);
    dphi[0] = -1.0;
    dphi[1] = 1.0;
};