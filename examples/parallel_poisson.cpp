#include <iostream>
#include "utilities/data_structures.hpp"
#include "utilities/export.hpp"
#include "utilities/norm.hpp"
#include "utilities/utils.hpp"
#include "utilities/mpi_util.hpp"
#include "quadrature/quadrature.hpp"
#include "mesh/mesh.hpp"
#include "fe/fem.hpp"
#include "solve/cg.hpp"


#include <set>
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

        MPI_Comm mpi_communicator;

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
        void exchange_shared();
        const size_t solve();
        void compute_errors();
        void output_results() const;

        std::vector<size_t> my_global_dofs;
        std::map<int, std::set<int>> shared_dofs;

        data_structures::parallel::SparseMatrix system_matrix;
        data_structures::parallel::Vector       system_rhs;
        //data_structures::parallel::Vector       solution;
        std::vector<double>                     solution;

        
    };

    Poisson1D::Poisson1D()
        :
        mpi_communicator(MPI_COMM_WORLD),
        n_mpi_processes(mpi_util::n_mpi_processes(mpi_communicator)),
        this_mpi_process(mpi_util::this_mpi_process(mpi_communicator)),
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
        mesh.partition(n_mpi_processes);    // split cells over processes
        mesh.distribute_dofs();             // distribute dofs over processes and mark shared dofs


        const std::vector<std::vector<size_t>> dof_distribution = mesh.get_distribution();
        my_global_dofs = dof_distribution[this_mpi_process];
        shared_dofs    = mesh.get_shared_dofs();

        // std::cout << "Rank: " << this_mpi_process << " Size: " << n_mpi_processes << "\n";

        // // print out map
        // for (const auto &pair : shared_dofs) 
        // {
        //     std::cout << pair.first << ": ";
        //     for (const auto &s : pair.second) 
        //     {
        //         std::cout << s << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << std::endl;

        // getchar();

        // for (auto cell = mesh.cell_begin(); cell != mesh.cell_end(); ++cell) {
        //     int subdomain = cell->get_subdomain();
        //     for (size_t i = 0; i < dofs_per_cell; ++i) {
        //         int dof = cell->vertex(i).global_index();
        //         shared_dofs[dof].insert(subdomain);
        //     }
        // }

        //system_rhs.assign(my_global_dofs.size(), 0.0);
        solution.assign(my_global_dofs.size(), 0.0);
        //solution.clear();
        system_rhs.clear();
        system_matrix.clear();

        // an std::map does not allocate memory sequentially, but I thought maybe 
        // if I initialize it like this, its placement in memory will be closer (this is probably stupid)
        for (const auto &my_dof : my_global_dofs)
        {
            //solution.insert({my_dof, 0.0});
            system_rhs.insert({my_dof, 0.0});
            system_matrix.insert(my_dof, my_dof, 0.0);
        }
        
        
        // const std::size_t n_dofs = mesh.n_cells() + 1;
        // system_matrix.resize(n_dofs, n_dofs, 0.0);
        // system_rhs.resize(n_dofs, 0.0);
        // solution.resize(n_dofs, 0.0);
    }

    void Poisson1D::assemble_system()
    {

        // if (!system_matrix.empty())
        //     throw std::runtime_error("Matrix is not empty");

        
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

    // void Poisson1D::exchange_shared()
    // {

    //     // Use a red-black coloring scheme to exchange shared dofs between processes
    //     // note: this algorithm assumes that maximally two processes share one dofs

    //     bool is_red = this_mpi_process % 2 == 0;

    //     for (const auto &dof : shared_dofs)
    //     {
    //         if (dof.second.size() > 1)
    //         {
    //             std::vector<double> send_buffer_rhs, send_buffer_matrix, receive_buffer_rhs, receive_buffer_matrix;

    //             const bool i_share_this_dof = dof.second.find(this_mpi_process) != dof.second.end();
                
    //             if (!i_share_this_dof)
    //                 continue;
                
    //             send_buffer_rhs.push_back(dof.first);
    //             send_buffer_rhs.push_back(system_rhs.at(dof.first));
    //             send_buffer_matrix.push_back(dof.first);
    //             send_buffer_matrix.push_back(system_matrix.at({dof.first, dof.first}));

    //             receive_buffer_rhs.resize(send_buffer_rhs.size());
    //             receive_buffer_matrix.resize(send_buffer_matrix.size());

    //             int other_process = -1;
    //             if (this_mpi_process == *dof.second.begin())
    //                 other_process = *std::prev(dof.second.end());
    //             else
    //                 other_process = *dof.second.begin();

    //             if (is_red) {

    //                 assert(other_process % 2 != 0);

    //                 // Red processes send first
    //                 MPI_Send(send_buffer_rhs.data(), send_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
    //                 MPI_Send(send_buffer_matrix.data(), send_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
    //                 MPI_Recv(receive_buffer_rhs.data(), receive_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);
    //                 MPI_Recv(receive_buffer_matrix.data(), receive_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);

    //             } else {
    //                 assert(other_process % 2 == 0);

    //                 // Black processes receive first
    //                 MPI_Recv(receive_buffer_rhs.data(), receive_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);
    //                 MPI_Recv(receive_buffer_matrix.data(), receive_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);
    //                 MPI_Send(send_buffer_rhs.data(), send_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
    //                 MPI_Send(send_buffer_matrix.data(), send_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
                    
    //             }

    //             for (size_t i = 0; i < receive_buffer_rhs.size(); i += 2)
    //             {
    //                 int dof = (int)receive_buffer_rhs[i];
    //                 double rhs_value = receive_buffer_rhs[i + 1];
    //                 int dof_matrix = (int)receive_buffer_matrix[i];
    //                 double matrix_value = receive_buffer_matrix[i + 1];
                    
    //                 system_rhs.at(dof) += rhs_value;
    //                 system_matrix.add(dof_matrix, dof_matrix, matrix_value);
    //             } 

    //             if (other_process == -1)
    //                 throw std::runtime_error("Other process not found");


    //         }
    //     }
    // }

    void Poisson1D::exchange_shared()
    {
        // Use a red-black coloring scheme to exchange shared dofs between processes
        // note: this algorithm assumes that maximally two processes share one dofs

        bool is_red = this_mpi_process % 2 == 0;

        for (const auto &dof : shared_dofs)
        {
            if (dof.second.size() > 1)
            {
                std::vector<double> send_buffer_rhs, send_buffer_matrix, receive_buffer_rhs, receive_buffer_matrix;

                const bool i_share_this_dof = dof.second.find(this_mpi_process) != dof.second.end();
                
                if (!i_share_this_dof)
                    continue;
                
                send_buffer_rhs.push_back(dof.first);
                send_buffer_rhs.push_back(system_rhs.at(dof.first));
                send_buffer_matrix.push_back(dof.first);
                send_buffer_matrix.push_back(system_matrix.at({dof.first, dof.first}));

                receive_buffer_rhs.resize(send_buffer_rhs.size());
                receive_buffer_matrix.resize(send_buffer_matrix.size());

                int other_process = -1;
                if (this_mpi_process == *dof.second.begin())
                    other_process = *std::prev(dof.second.end());
                else
                    other_process = *dof.second.begin();

                if (is_red) {

                    assert(other_process % 2 != 0);

                    // Red processes send first
                    MPI_Send(send_buffer_rhs.data(), send_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
                    MPI_Send(send_buffer_matrix.data(), send_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
                    MPI_Recv(receive_buffer_rhs.data(), receive_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);
                    MPI_Recv(receive_buffer_matrix.data(), receive_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);

                } else {
                    assert(other_process % 2 == 0);

                    // Black processes receive first
                    MPI_Recv(receive_buffer_rhs.data(), receive_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);
                    MPI_Recv(receive_buffer_matrix.data(), receive_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator, MPI_STATUS_IGNORE);
                    MPI_Send(send_buffer_rhs.data(), send_buffer_rhs.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);
                    MPI_Send(send_buffer_matrix.data(), send_buffer_matrix.size(), MPI_DOUBLE, other_process, 0, mpi_communicator);

                }

                for (size_t i = 0; i < receive_buffer_rhs.size(); i += 2)
                {
                    int dof = (int)receive_buffer_rhs[i];
                    double rhs_value = receive_buffer_rhs[i + 1];
                    int dof_matrix = (int)receive_buffer_matrix[i];
                    double matrix_value = receive_buffer_matrix[i + 1];
                    
                    if (is_red)
                    {
                        system_rhs.at(dof) += rhs_value;
                        system_matrix.add(dof_matrix, dof_matrix, matrix_value);
                    }
                    else
                    {
                        system_rhs.erase(dof);
                        system_matrix.erase({dof_matrix, dof_matrix});
                        //solution.erase(solution.end()); // resize the solution vector
                        solution.resize(solution.size() - 1);
                    }

                } 

                if (other_process == -1)
                    throw std::runtime_error("Other process not found");
            }
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
        exchange_shared();      // exchange shared dofs between processes


        const size_t cg_iterations = solve();
        std::cout << "CG iterations: " << cg_iterations << std::endl;

        if (this_mpi_process == 0)
        {
            for (const auto &s : solution)
                std::cout << s << " ";
            std::cout << std::endl;
            
        }

        // int local_size = my_global_dofs.size();

        // // Gather the sizes of all local solution vectors
        // std::vector<int> local_sizes(n_mpi_processes);
        // MPI_Allgather(&local_size, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        // // Calculate the total size of the global solution vector and the displacements for each local solution vector
        // int global_size = 0;
        // std::vector<int> displacements(n_mpi_processes);
        // for (int i = 0; i < n_mpi_processes; ++i) {
        //     displacements[i] = global_size;
        //     global_size += local_sizes[i];
        // }

        // // Create a vector to hold the global solution
        // std::vector<double> global_solution(global_size);

        // // Gather all local solution vectors into the global solution vector
        // MPI_Gatherv(solution.data(), local_size, MPI_DOUBLE, global_solution.data(), local_sizes.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);


        // if (this_mpi_process == 0)
        // {
        //     for (const auto &s : global_solution)
        //         std::cout << s << " ";
        //     std::cout << std::endl;

        // }
            
    }




}   // namespace parallel_poisson




int main(int argc, char **argv)
{

    using namespace parallel_poisson;

    //MPI::Environment env(argc, argv);
    mpi_util::MPIUtil mpi_env(argc, argv);
    //mpi_env.init(argc, argv);

    Poisson1D poisson;

    //std::cout << "Rank: " << mpi_env.get_rank() << " Size: " << mpi_env.get_size() << "\n";
    

    poisson.run();

    //poisson.mesh_info();

    //poisson.run();


    return 0;
}
