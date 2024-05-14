#include <iostream>
#include "utilities/data_structures.hpp"
#include "utilities/export.hpp"
#include "utilities/norm.hpp"
#include "utilities/utils.hpp"
#include "quadrature/quadrature.hpp"
#include "mesh/mesh.hpp"
#include "fe/fem.hpp"
#include "solve/cg.hpp"

#include <chrono>



namespace parallel_poisson 
{

    class Poisson1D
    {
    public:
        
        Poisson1D();
        
        void run();
        void mesh_info(const bool detailed = false);

    private:

        static const double f(const Point<1> &x);
        static const double u(const Point<1> &x);

        //MPI_Comm mpi_communicator;

        const std::size_t n_mpi_processes;
        const std::size_t this_mpi_process;

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


        void setup_system();
        void assemble_system();
        const size_t solve();
        void compute_errors();
        void output_results() const;

        std::vector<size_t> my_global_dofs;

        data_structures::parallel::SparseMatrix system_matrix;
        data_structures::parallel::Vector       system_rhs;
        data_structures::parallel::Vector       solution;

        
    };

    Poisson1D::Poisson1D()
        :
        // mpi_communicator(MPI_COMM_WORLD),
        // n_mpi_processes(MPI::n_mpi_processes(mpi_communicator)),
        // this_mpi_process(MPI::this_mpi_process(mpi_communicator)),
        n_mpi_processes(3),
        this_mpi_process(0),
        mesh(0, 1., 10),
        psi(),
        qr(quadrature::midpoint),
        dofs_per_cell(psi.ndof),
        n_quad_pts_per_cell(qr.n)
    {
        bc.g = u;
        bc.lbs = {1, 2};
        bc.set_dirichlet = true;
    }

    void Poisson1D::mesh_info(const bool detailed) 
    {
        if (this_mpi_process == 0) 
        {
            mesh.mesh_info(detailed);
        }
    
    }

    // Right hand side function
    const double Poisson1D::f(const Point<1> & x) 
    {
        return (M_PI*M_PI*(361*cos((19*M_PI*x[0])/10) - 441*cos((21*M_PI*x[0])/10)))/100;
    }
    // Exact solution
    const double Poisson1D::u(const Point<1> & x) 
    {
        return 2*sin(2*M_PI*x[0])*sin(M_PI*x[0]/10) + 10;
    }

    void Poisson1D::compute_stiffness_on_cell(
        const typename Mesh1D::cell_iterator &cell,
        data_structures::serial::DenseMatrix &Ak) 
    {

        Ak.reinit(0.);

        // Holder for evaluations of the gradient of psi
        data_structures::serial::DenseMatrix dpsi_vals(dofs_per_cell, dim);
        
        const double measure = cell->get_measure();

        // Loop over trial function dofs
        for (int i = 0; i < dofs_per_cell; ++i) 
        {        
            // Loop over test function dofs
            for (int j = 0; j < dofs_per_cell; ++j) 
            {    
                // Integrate shape functions of current dofs over the cell 
                for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) 
                {
                    const Point<1> xq(qr[ipq].node);   // quadrature point in reference element
                    psi.eval_d(*cell, xq, dpsi_vals);    // evaluate gradient of psi at xq -> store in dpsi_vals

                    for (int dm = 0; dm < dim; dm++)     // loop over space dimensions
                    { 
                        Ak(i, j) += qr[ipq].weight * measure * dpsi_vals(i, dm) * dpsi_vals(j, dm);
                    }       
                }
            }
        }
    }

    void Poisson1D::compute_rhs_on_cell(
        const typename Mesh1D::cell_iterator &cell,
        Vector &fk)
    {

        fk.reinit(0.0);

        Vector psi_vals(dofs_per_cell);    // container for evaluations of psi

        // Get the measure of the element
        const double measure = cell->get_measure();

        // Loop over dofs
        for (int i = 0; i < dofs_per_cell; i++) 
        {
        
            // Loop over quadrature points
            for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) 
            {

                Point<dim> x;
                const Point<dim> xq(qr[ipq].node); // quadrature point in reference element
                cell->map_to_physical(xq, x);      // map xq to quadrature point in physical element x
                psi.eval(xq, psi_vals);            // evaluate psi at xq -> store in psi_vals

                fk[i] += qr[ipq].weight * measure * f(x) * psi_vals[i];
            }

        }

    }

    void Poisson1D::get_boundary_data(
        const typename Mesh1D::cell_iterator &cell,
        const std::vector<size_t> &loc2glb,
        std::map<int, double> &boundary_data)
    {
        for (int i = 0; i < loc2glb.size(); i++) 
        {
            if (bc.set_dirichlet) 
            {
                // Add data if dof is a boundary dof
                if (std::find(bc.lbs.begin(), bc.lbs.end(), cell->vertex(i).boundary_label()) != bc.lbs.end())
                    boundary_data.emplace(cell->vertex(i).global_index(), bc.g(cell->vertex(i)));       
            }
        }

    }

    void Poisson1D::setup_system()
    {
        mesh.partition(n_mpi_processes);

        const std::vector<std::vector<size_t>> dof_distribution
        = mesh.get_distribution();  // index_distribution[p][i] = global index of local dof i on process p

        my_global_dofs = dof_distribution[this_mpi_process];

        system_rhs.assign(my_global_dofs.size(), 0.0);
        solution.assign(my_global_dofs.size(), 0.0);
        system_matrix.clear();

        
        
        // const std::size_t n_dofs = mesh.n_cells() + 1;
        // system_matrix.resize(n_dofs, n_dofs, 0.0);
        // system_rhs.resize(n_dofs, 0.0);
        // solution.resize(n_dofs, 0.0);
    }

    void Poisson1D::assemble_system()
    {

        if (!system_matrix.empty())
            throw std::runtime_error("Matrix is not empty");

        
        data_structures::serial::DenseMatrix Ak(dofs_per_cell, dofs_per_cell);    
        data_structures::serial::Vector fk(dofs_per_cell);
        std::map<int, double> boundary_data;

        // Loop over all cells in the mesh
        for (auto cell = mesh.cell_begin(); cell != mesh.cell_end(); ++cell) 
            if (cell->get_subdomain() == this_mpi_process) 
            {
                std::vector<size_t> loc2glb(dofs_per_cell);
                for (size_t i = 0; i < dofs_per_cell; ++i) 
                    loc2glb[i] = cell->vertex(i).global_index();

                compute_stiffness_on_cell(cell, Ak);      
                compute_rhs_on_cell(cell, fk);
                get_boundary_data(cell, loc2glb, boundary_data);    
                utilities::parallel_dofs::distribute_local_to_global(
                    Ak, 
                    fk,
                    system_matrix, 
                    system_rhs,
                    loc2glb,
                    boundary_data);

            }

    }

    const size_t Poisson1D::solve()
    {
        const size_t max_iter = 1000;
        const double tol = 1e-10;
        return solve::parallel::cg(system_matrix, system_rhs, solution, max_iter, tol);    
    }

    void Poisson1D::run()
    {
        setup_system();
        assemble_system();
    }




}   // namespace parallel_poisson




int main(int argc, char **argv)
{
    //MPI::Environment env(argc, argv);

    parallel_poisson::Poisson1D poisson;

    poisson.run();

    poisson.mesh_info();

    //poisson.run();

    return 0;
}
