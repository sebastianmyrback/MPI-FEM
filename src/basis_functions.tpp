#include "basis_functions.hpp"
template<int degree>
void lagrange_1d<degree>::eval(const Rn &x, std::vector<double> &phi) const {
    // x - quadrature point in reference element 

    phi.resize(this->ndof);
    phi[0] = 1.0 - x[0];
    phi[1] = x[0];
};

template<int degree>
void lagrange_1d<degree>::eval_d(const FE &K, const Rn &x, std::vector<std::vector<double>> &dphi) const {
    // K - physical element
    // x - quadrature point in reference element
    // dphi - matrix of size ndof x D holding the derivatives of the basis functions

    assert(dphi.size() == this->ndof);
    assert(dphi[0].size() == 1);

    //const std::vector<double> vertices = {K(0).x[0], K(1).x[0]};
    const std::vector<double> vertices = {K.elem_vertices[0]->x[0], K.elem_vertices[1]->x[0]};

    for (int i = 0; i < this->ndof; i++) {
        int i1 = (i+1)%2;
        dphi[i][0] = 1. / (vertices[i] - vertices[i1]);
    }

};