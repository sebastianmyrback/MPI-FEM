#include <iostream>
#include "mesh.hpp"
#include <assert.h>
#include "fe.hpp"
#include "problem.hpp"
#include "lagrange.hpp"

int main() {
    
    // Create a mesh object
    const int n = 10;     // number of elements
    const int n_threads = 1;
    const double a = -2.0, b = 2.0;

    mesh1d Th(a, b, n);

    QuadratureRule<1> trapezoidal(2, {QuadraturePoint<1>(Rd<1>({0.0}), 0.5), QuadraturePoint<1>(Rd<1>({1.0}), 0.5)});

    lagrange_1d<1> psi;

    problem<mesh1d> prob(n_threads);

    prob.assemble_FEM_matrix(Th, trapezoidal, psi, 0., 1.);


    //gnuplot::save(Th);

    return 0;
}