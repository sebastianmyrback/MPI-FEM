#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <map>
#include <functional>
#include "basis_functions.hpp"
#include "quadrature.hpp"

template <int d>
struct dirichlet_bc {

    typedef Rd<d> Rn;

    std::function<double(const Rn &)> g;    // Dirichlet boundary function
    std::vector<int> lbs;                   // labels of Dirichlet boundaries
    bool set_dirichlet = false;             // if false, no Dirichlet boundary condition is set strongly


};


template <typename mesh>
class problem {

    typedef std::map<std::pair<int, int>, double> matrix;
    typedef typename mesh::Rn Rn;
    typedef typename mesh::Element Quad;

    std::shared_ptr<mesh> Th;

public:

    matrix mat;              // System matrix 

    std::vector<double> rhs; // Right hand side vector

    const int thread_count;  // Number of threads to use

    const int n_dofs;        // Total number of degrees of freedom of the problem

    problem(const int thread_count, std::shared_ptr<mesh> Th)
        : thread_count(thread_count), n_dofs(0), Th(std::move(Th)) {}

    template <int d, int degree>
    void assemble_FEM_matrix(const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const double mass, const double stiffness, const dirichlet_bc<d> & bc);
    
    template <int d, int degree>
    void assemble_FEM_matrix(const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const double mass, const double stiffness) {
        dirichlet_bc<d> bc;    // set_dirichlet defaults to false
        assemble_FEM_matrix(qr, psi, mass, stiffness, bc);
    };

    template <int d, int degree>
    void assemble_rhs(const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const std::function<double(const typename mesh::Rn &)> & f);

    void set_dirichlet(const std::function<double(const typename mesh::Rn &)> & g);

};



#include "problem.tpp"

#endif