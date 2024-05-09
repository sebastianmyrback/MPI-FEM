#pragma once


#include <map>
#include <unordered_map>
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

    std::shared_ptr<mesh> Th;

public:

    Matrix mat;              // System matrix 
    Vector rhs;              // Right hand side vector

    const int thread_count;  // Number of threads to use
    const int n_dofs;        // Total number of degrees of freedom of the problem

    FEM(const int thread_count, std::shared_ptr<mesh> Th)
        : thread_count(thread_count), n_dofs(0), Th(std::move(Th)) {}

    // Assembly methods

    template <int dim, int degree>
    void compute_stiffness_on_cell(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        const QuadratureRule<dim> &qr,
        const BasisFunction<dim, degree> &psi,
        std::vector<std::vector<double>> &Ak); 

    template <int dim, int degree>
    void compute_rhs_on_cell(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        const QuadratureRule<dim> &qr,
        const BasisFunction<dim, degree> &psi,
        const double f(const Point<dim> &),
        std::vector<double> &fk);

    template <int dim>
    void get_boundary_data(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        const DirichletBC<dim> &bc,
        std::unordered_map<int, double> &boundary_data);

    void distribute_local_to_global(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        std::vector<std::vector<double>> &Ak,
        std::vector<double> &fk,
        const std::unordered_map<int, double> &boundary_data);

    template <int dim, int degree>
    void assemble_stiffness_system(
        const QuadratureRule<dim> &qr, 
        const BasisFunction<dim, degree> & psi, 
        const double f(const Point<dim> &),
        const DirichletBC<dim> & bc);

    template <int dim, int degree>
    void assemble_stiffness_system(
        const QuadratureRule<dim> &qr, 
        const BasisFunction<dim, degree> & psi, 
        const double f(const Point<dim> &))
    {
        DirichletBC<dim> bc;    // set_dirichlet defaults to false
        assemble_stiffness_system(qr, psi, f, bc);
    };

    
};



#include "fem.tpp"

