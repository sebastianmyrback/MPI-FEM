#ifndef FEM_HPP
#define FEM_HPP

#include <map>
#include <functional>
#include "basis_functions.hpp"
#include "quadrature.hpp"

template <int d>
struct DirichletBC {

    typedef Point<d> Rn;

    std::function<double(const Rn &)> g;    // Dirichlet boundary function
    std::vector<int> lbs;                   // labels of Dirichlet boundaries
    bool set_dirichlet = false;             // if false, no Dirichlet boundary condition is set strongly


};


template <typename mesh>
class FEM {

    typedef std::map<std::pair<int, int>, double> Matrix;
    typedef std::vector<double> Vector;
    typedef typename mesh::Rn Rn;

    std::shared_ptr<mesh> Th;

public:

    Matrix mat;              // System matrix 
    Vector rhs;              // Right hand side vector

    const int thread_count;  // Number of threads to use
    const int n_dofs;        // Total number of degrees of freedom of the problem

    FEM(const int thread_count, std::shared_ptr<mesh> Th)
        : thread_count(thread_count), n_dofs(0), Th(std::move(Th)) {}

    template <int dim, int degree>
    void assemble_FEM_matrix(const QuadratureRule<dim> &qr, const BasisFunction<dim, degree> & psi, const double mass, const double stiffness, const DirichletBC<dim> & bc);
    
    template <int dim, int degree>
    void assemble_FEM_matrix(const QuadratureRule<dim> &qr, const BasisFunction<dim, degree> & psi, const double mass, const double stiffness) {
        DirichletBC<dim> bc;    // set_dirichlet defaults to false
        assemble_FEM_matrix(qr, psi, mass, stiffness, bc);
    };

    template <int d, int degree>
    void assemble_rhs(const QuadratureRule<d> &qr, const BasisFunction<d, degree> & psi, const std::function<double(const typename mesh::Rn &)> & f);

    void set_dirichlet(const std::function<double(const typename mesh::Rn &)> & g);

};



#include "fem.tpp"

#endif  // FEM_HPP