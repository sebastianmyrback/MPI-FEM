// Write assembly routines for the stiffness matrix

// Path: src/fem/Assemble.hpp

#ifndef ASSEMBLE_HPP
#define ASSEMBLE_HPP

#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>

#include "mesh.hpp"
#include "CG.hpp"
#include "Lagrange.hpp"

// Assemble the stiffness matrix
// A is a square matrix of size n x n
// A is stored as a 1D array of size n*n
// A[i*n + j] is the element in the ith row and jth column

// The function returns the number of iterations
// The function modifies the matrix A

template <typename mesh_t>
int assemble_stiffness_matrix(double *A, mesh_t &Th) {
    
    const int d = Th.D;
    typedef finite_element<d> FiniteElement;
    
    // Get the number of vertices
    int n = Th.nv;

    // Get the number of elements
    int ne = Th.nk;

    // Get the degree of the finite element
    int deg = P1Lagrange1D::deg;

    // Get the number of degrees of freedom
    int ndof = P1Lagrange1D::ndof;

    // Initialize the stiffness matrix
    for (int i = 0; i < n*n; i++) {
        A[i] = 0.0;
    }

    // Loop over the elements
    for (int k = 0; k < ne; k++) {
        // Get the element
        element<1> &elem = Th[k];

        // Get the vertices of the element
        vertex &v0 = elem(0);
        vertex &v1 = elem(1);

        // Get the coordinates of the vertices
        double x0 = v0.x;
        double x1 = v1.x;

        // Get the length of the element
        double h = x1 - x0;

        // Get the quadrature points and weights
        QuadratureFormula<d> quad;

        std::vector<double> q = {0.5};
        std::vector<double> w = {1.0};

        // Loop over the quadrature points
        for (int q = 0; q < 1; q++) {
            // Get the quadrature point
            double x = 0.5*(x0 + x1) + 0.5*h*q;

            // Evaluate the basis functions
            std::vector<double> phi;
            P1Lagrange1D::eval(x, phi);

            // Evaluate the derivatives of the basis functions
            std::vector<double> dphi;
            P1Lagrange1D::eval_d(x, dphi);

            // Loop over the degrees of freedom
            for (int i = 0; i < ndof; i++) {
                for (int j = 0; j < ndof; j++) {
                    A[elem(i)*n + elem(j)] += h*w[q]*dphi[i]*dphi[j];
                }
            }
        }
    }

    return 0;
}


#endif