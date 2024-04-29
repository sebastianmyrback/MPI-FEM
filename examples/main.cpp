#include <iostream>
#include "mesh.hpp"
#include <assert.h>
#include "fe.hpp"
#include "assemble.hpp"
#include "lagrange.hpp"

int main() {
    
    // Create a mesh object
    int n = 10;     // number of elements
    double a = 0.0, b = 1.0;

    mesh1d Th(a, b, n);

    std::vector<double> x = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    p1_lagrange_1d Vh;

    problem<mesh1d> prob(1);

    prob.assemble_FEM_matrix(Th, Vh, 0., 1.);




    //gnuplot::save(Th);

    return 0;
}