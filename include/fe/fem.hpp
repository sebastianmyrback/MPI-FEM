#pragma once

#include <map>
#include <functional>
#include "../fe/basis_functions.hpp"
#include "../quadrature/quadrature.hpp"
#include "../utilities/data_structures.hpp"
#include "../utilities/utils.hpp"


using namespace data_structures::serial;

template <typename mesh>
class FEM {

    const mesh* Th;

public:

    SparseMatrix mat;              // System matrix 
    Vector rhs;              // Right hand side vector

    const int thread_count;  // Number of threads to use
    const int n_dofs;        // Total number of degrees of freedom of the problem

    FEM(const int thread_count, const mesh *_Th)
        : thread_count(thread_count), n_dofs(0), Th(std::move(_Th)) {}

    // Assembly methods

    // Compute the local stiffness matrix on a cell
    template <int dim, int degree>
    void compute_stiffness_on_cell(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        const quadrature::QuadratureRule<dim> &qr,
        const BasisFunction<dim, degree> &psi,
        DenseMatrix &Ak); 

    // Compute the local rhs vector on a cell
    template <int dim, int degree>
    void compute_rhs_on_cell(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        const quadrature::QuadratureRule<dim> &qr,
        const BasisFunction<dim, degree> &psi,
        const double f(const Point<dim> &),
        Vector &fk);

    // This function checks if the current cell
    // contains any boundary dofs and adds the corresponding
    // global dof indices and boundary values to the boundary_data map
    template <int dim>
    void get_boundary_data(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        const utilities::DirichletBC<dim> &bc,
        std::map<int, double> &boundary_data);

    // This function adds the local stiffness matrix Ak and local rhs vector
    // fk into the global matrix mat and global rhs vector rhs.
    // Moreover, it applies strong Dirichlet boundary conditions to the system
    // if the user has specified so in the DirichletBC object
    void distribute_local_to_global(
        const typename mesh::cell_iterator &cell,
        const std::vector<int> &loc2glb,
        DenseMatrix &Ak,
        Vector &fk,
        const std::map<int, double> &boundary_data);


    // The following function assembles the system using the functions above
    template <int dim, int degree>
    void assemble_stiffness_system(
        const quadrature::QuadratureRule<dim> &qr, 
        const BasisFunction<dim, degree> & psi, 
        const double f(const Point<dim> &),
        const utilities::DirichletBC<dim> & bc);

    template <int dim, int degree>
    void assemble_stiffness_system(
        const quadrature::QuadratureRule<dim> &qr, 
        const BasisFunction<dim, degree> & psi, 
        const double f(const Point<dim> &))
    {
        utilities::DirichletBC<dim> bc;    // set_dirichlet defaults to false
        assemble_stiffness_system(qr, psi, f, bc);
    };

    
};



#include "fem.tpp"

