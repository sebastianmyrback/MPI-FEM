#pragma once

#include <iostream>
#include <set>
#include "utilities/data_structures.hpp"
#include "utilities/export.hpp"
#include "utilities/norm.hpp"
#include "utilities/utils.hpp"
#include "utilities/mpi_util.hpp"
#include "quadrature/quadrature.hpp"
#include "mesh/mesh.hpp"
#include "fe/fem.hpp"
#include "solve/cg.hpp"




namespace parallel_poisson 
{

    class Poisson1D
    {
    public:
        
        Poisson1D(const double a, const double b, const int nintervals);
        
        void run();
        void setup_system();
        void assemble_system();
        void exchange_shared();
        const size_t solve();
        const double L2H1norm(const double l2, const double h1);    // not implemented
        void output_solution(const std::string &filename) const;

        void mesh_info(const bool detailed = false);

    private:

        static const double f(const Point<1> &x);   // rhs function
        static const double u(const Point<1> &x);   // exact solution

        MPI_Comm mpi_communicator;
        const size_t n_mpi_processes;
        const size_t this_mpi_process;

        Mesh1D                              mesh;
        const quadrature::QuadratureRule<1> qr;
        const P1Lagrange1D<1>               psi;
        utilities::DirichletBC<1>           bc;

        const int dofs_per_cell;
        const int n_quad_pts_per_cell;
        constexpr static int dim = Mesh1D::dim;


        void compute_stiffness_on_cell(
            const typename Mesh1D::cell_iterator &cell, 
            data_structures::serial::DenseMatrix &Ak);

        void compute_rhs_on_cell(
            const typename Mesh1D::cell_iterator &cell, 
            data_structures::serial::Vector &fk);

        void get_boundary_data(
            const typename Mesh1D::cell_iterator &cell, 
            const std::vector<size_t> &loc2glb,
            std::map<int, double> &boundary_data);   

        std::vector<size_t> my_global_dofs;
        std::map<int, std::set<int>> shared_dofs;

        data_structures::parallel::SparseMatrix system_matrix;
        data_structures::parallel::Vector       system_rhs;
        std::vector<double>                     solution;

        
    };

}   // namespace parallel_poisson



namespace serial_poisson 
{

    class Poisson1D
    {
    public:
        
        Poisson1D(const double a, const double b, const int nintervals);
        
        void run();
        void setup_system();
        void assemble_system();
        const size_t solve();
        const double L2H1norm(const double l2, const double h1);    // not implemented
        void output_solution(const std::string &filename) const;

        void mesh_info(const bool detailed = false);

    private:

        static const double f(const Point<1> &x);   // rhs function
        static const double u(const Point<1> &x);   // exact solution

        Mesh1D                              mesh;
        const quadrature::QuadratureRule<1> qr;
        const P1Lagrange1D<1>               psi;
        utilities::DirichletBC<1>           bc;

        const int dofs_per_cell;
        const int n_quad_pts_per_cell;
        constexpr static int dim = Mesh1D::dim;


        void compute_stiffness_on_cell(
            const typename Mesh1D::cell_iterator &cell, 
            data_structures::serial::DenseMatrix &Ak);

        void compute_rhs_on_cell(
            const typename Mesh1D::cell_iterator &cell, 
            data_structures::serial::Vector &fk);

        void get_boundary_data(
            const typename Mesh1D::cell_iterator &cell, 
            const std::vector<size_t> &loc2glb,
            std::map<int, double> &boundary_data);   



        std::vector<size_t> my_global_dofs;
        std::map<int, std::set<int>> shared_dofs;

        data_structures::serial::SparseMatrix system_matrix;
        data_structures::serial::Vector       system_rhs;
        std::vector<double>                     solution;

        
    };

}   // namespace serial_poisson

