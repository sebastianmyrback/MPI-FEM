#ifndef PROBLEM_HPP
#define PROBLEM_HPP

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

    template <int d, int degree>
    void assemble_FEM_matrix(const mesh & Th, const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const double alpha, const double beta);

};



#include "../src/problem.tpp"

#endif