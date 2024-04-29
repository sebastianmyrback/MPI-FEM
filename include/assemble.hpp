// Write assembly routines for the stiffness matrix

// Path: src/fem/Assemble.hpp

#ifndef ASSEMBLE_HPP
#define ASSEMBLE_HPP

#include <map>
#include "mesh.hpp"
#include "lagrange.hpp"
#include "quadrature.hpp"

template <typename mesh>
class problem {

    typedef std::map<std::pair<int, int>, double> matrix;

public:

    matrix mat;              // System matrix 

    std::vector<double> rhs; // Right hand side vector

    const int thread_count;  // Number of threads to use

    const int n_dofs;        // Number of degrees of freedom

    problem(const int thread_count) : thread_count(thread_count), n_dofs(0) {}

    void assemble_FEM_matrix(const mesh & Th, const p1_lagrange_1d & Vh, const double alpha, const double beta);

};



#include "../src/assemble.tpp"

#endif