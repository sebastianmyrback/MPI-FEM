#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <map>
#include <functional>
#include "basis_functions.hpp"
#include "quadrature.hpp"


template <typename mesh>
class problem {

    typedef std::map<std::pair<int, int>, double> matrix;
    typedef typename mesh::Rn Rn;
    typedef typename mesh::elem elem;

public:

    matrix mat;              // System matrix 

    std::vector<double> rhs; // Right hand side vector

    const int thread_count;  // Number of threads to use

    const int n_dofs;        // Number of degrees of freedom

    problem(const int thread_count) : thread_count(thread_count), n_dofs(0) {}

    template <int d, int degree>
    void assemble_FEM_matrix(const mesh & Th, const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const double mass, const double stiffness);

    template <int d, int degree>
    void assemble_rhs(const mesh & Th, const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const std::function<double(const typename mesh::Rn &)> & f);

    void set_dirichlet(const mesh & Th, const std::function<double(const typename mesh::Rn &)> & g);

};



#include "../src/problem.tpp"

#endif