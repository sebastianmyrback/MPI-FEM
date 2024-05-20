#include "../include/poisson.hpp"


namespace serial_poisson 
{
    Poisson1D::Poisson1D(const double a, const double b, const int nintervals)
        : mesh(a, b, nintervals),
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
        mesh.mesh_info(detailed);
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
        std::vector<std::vector<double>> dpsi_vals(dofs_per_cell, std::vector<double>(dim, 0.));
        
        const double measure = cell->get_measure();

        // Integrate shape functions of current dofs over the cell 
        for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) 
        {
            const Point<1> xq(qr[ipq].node);   // quadrature point in reference element
            psi.eval_d(*cell, xq, dpsi_vals);    // evaluate gradient of psi at xq -> store in dpsi_vals
            
            // Loop over trial function dofs
            for (int i = 0; i < dofs_per_cell; ++i) 
            {        
                // Loop over test function dofs
                for (int j = 0; j < dofs_per_cell; ++j) 
                {    
                    for (int dm = 0; dm < dim; dm++)     // loop over space dimensions
                    { 
                        Ak(i, j) += qr[ipq].weight * measure * dpsi_vals[i][dm] * dpsi_vals[j][dm];
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

        //Vector psi_vals(dofs_per_cell);    // container for evaluations of psi
        std::vector<double> psi_vals(dofs_per_cell);    // container for evaluations of psi

        // Get the measure of the element
        const double measure = cell->get_measure();

        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) 
        {
            Point<dim> x;
            const Point<dim> xq(qr[ipq].node); // quadrature point in reference element
            cell->map_to_physical(xq, x);      // map xq to quadrature point in physical element x
            psi.eval(xq, psi_vals);            // evaluate psi at xq -> store in psi_vals

            // Loop over dofs
            for (int i = 0; i < dofs_per_cell; i++) 
            {
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
        system_rhs.assign(mesh.get_nverts(), 0.0);
        solution.assign(mesh.get_nverts(), 0.0);
    }

    void Poisson1D::assemble_system()
    {

        data_structures::serial::DenseMatrix Ak(dofs_per_cell, dofs_per_cell);    
        data_structures::serial::Vector fk(dofs_per_cell);
        std::map<int, double> boundary_data;

        // Loop over all cells in the mesh
        for (auto cell = mesh.cell_begin(); cell != mesh.cell_end(); ++cell) 
        {
            std::vector<size_t> loc2glb(dofs_per_cell);
            for (size_t i = 0; i < dofs_per_cell; ++i) 
                loc2glb[i] = cell->vertex(i).global_index();

            compute_stiffness_on_cell(cell, Ak);      
            
            compute_rhs_on_cell(cell, fk);
            
            get_boundary_data(cell, loc2glb, boundary_data);    
            
            utilities::distribute_local_to_global(
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
        return solve::serial::cg(system_matrix, system_rhs, solution, max_iter, tol);    
    }

    const double Poisson1D::L2H1norm(const double l2, const double h1)
    {
        
        assert(false);  // not correctly implemented yet

    }

    void Poisson1D::output_solution(const std::string &filename) const
    {
        // Let each process write out "solution" to its own file
        std::ofstream file;
        file.open(filename + ".dat");
        for (size_t i = 0; i < solution.size(); i++)
            file << solution[i] << std::endl;
        file.close();

        std::cout << "Solution written to " << filename + ".dat" << std::endl;
    }
    
    void Poisson1D::run()
    {
        setup_system();
        assemble_system();
        const size_t cg_iterations = solve();
        std::cout << "CG iterations: " << cg_iterations << std::endl;
        output_solution("solution");

        // if (this_mpi_process == 0) 
        // {
            
        //     // print soluion
        //     for (size_t i = 0; i < solution.size(); i++)
        //         std::cout << solution[i] << std::endl;
        //     std::cout << std::endl;   
        
        // }

    }

}   // namespace serial_poisson



namespace parallel_poisson 
{
    Poisson1D::Poisson1D(const double a, const double b, const int nintervals)
        :
        mpi_communicator(MPI_COMM_WORLD),
        n_mpi_processes(mpi_util::n_mpi_processes(mpi_communicator)),
        this_mpi_process(mpi_util::this_mpi_process(mpi_communicator)),
        mesh(a, b, nintervals, n_mpi_processes),
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
        std::vector<std::vector<double>> dpsi_vals(dofs_per_cell, std::vector<double>(dim, 0.));
        
        const double measure = cell->get_measure();

        // Integrate shape functions of current dofs over the cell 
        for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) 
        {
            const Point<1> xq(qr[ipq].node);   // quadrature point in reference element
            psi.eval_d(*cell, xq, dpsi_vals);    // evaluate gradient of psi at xq -> store in dpsi_vals

            // Loop over trial function dofs
            for (int i = 0; i < dofs_per_cell; ++i) 
            {        
                // Loop over test function dofs
                for (int j = 0; j < dofs_per_cell; ++j) 
                {                
                    for (int dm = 0; dm < dim; dm++)     // loop over space dimensions
                    { 
                        Ak(i, j) += qr[ipq].weight * measure * dpsi_vals[i][dm] * dpsi_vals[j][dm];
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

        //Vector psi_vals(dofs_per_cell);    // container for evaluations of psi
        std::vector<double> psi_vals(dofs_per_cell);    // container for evaluations of psi

        // Get the measure of the element
        const double measure = cell->get_measure();

        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) 
        {
            Point<dim> x;
            const Point<dim> xq(qr[ipq].node); // quadrature point in reference element
            cell->map_to_physical(xq, x);      // map xq to quadrature point in physical element x
            psi.eval(xq, psi_vals);            // evaluate psi at xq -> store in psi_vals

            // Loop over dofs
            for (int i = 0; i < dofs_per_cell; i++) 
                fk[i] += qr[ipq].weight * measure * f(x) * psi_vals[i];
            

        }

    }

    void Poisson1D::get_boundary_data(
        const typename Mesh1D::cell_iterator &cell,
        const std::vector<size_t> &loc2glb,
        std::map<int, double> &boundary_data)
    {
        for (int i = 0; i < loc2glb.size(); i++) 
            if (bc.set_dirichlet)     
                if (std::find(bc.lbs.begin(), bc.lbs.end(), cell->vertex(i).boundary_label()) != bc.lbs.end())
                    boundary_data.emplace(cell->vertex(i).global_index(), bc.g(cell->vertex(i)));       
            
    }

    void Poisson1D::setup_system()
    {
        //mesh.partition(n_mpi_processes);    // split cells over processes
        //mesh.distribute_dofs();             // distribute dofs over processes and mark shared dofs

        const std::vector<std::vector<size_t>> dof_distribution = mesh.get_distribution();
        my_global_dofs = dof_distribution[this_mpi_process];
        shared_dofs    = mesh.get_shared_dofs();

        solution.assign(my_global_dofs.size(), 0.0);
        system_rhs.clear();
        system_matrix.clear();

        // an std::map does not allocate memory sequentially, but I thought maybe 
        // if I initialize it like this, its placement in memory will be closer (this is probably stupid)
        for (const auto &my_dof : my_global_dofs)
        {
            system_rhs.insert({my_dof, 0.0});
            system_matrix.insert(my_dof, my_dof, 0.0);
        }
        
    }

    void Poisson1D::assemble_system()
    {

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
                
                utilities::distribute_local_to_global(
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
    //             send_buffer_matrix.push_back(system_matrix.at({dof.first, dof.first}));
    //             send_buffer_matrix.push_back(system_matrix.at({dof.first, dof.first + 1 - 2 * (dof.first > my_global_dofs[0])}));    // send the off-diagonal entry 

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
    //                 double matrix_value = receive_buffer_matrix[i];
    //                 double matrix_value_next = receive_buffer_matrix[i + 1];
                    
    //                 if (is_red)
    //                 {
    //                     system_rhs.at(dof) += rhs_value;
    //                     system_matrix.add(dof, dof, matrix_value);
    //                     system_matrix.add(dof, dof + 1 - 2 * (dof <= my_global_dofs[0]), matrix_value_next);
    //                 }
    //                 else
    //                 {
    //                     system_rhs.erase(dof);
    //                     system_matrix.erase({dof, dof});
    //                     system_matrix.erase({dof, dof + 1 - 2 * (dof > my_global_dofs[0])});
    //                 }

    //             } 
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
                send_buffer_matrix.push_back(system_matrix.at({dof.first, dof.first}));
                send_buffer_matrix.push_back(system_matrix.at({dof.first, dof.first + 1 - 2 * (dof.first > my_global_dofs[0])}));    // send the off-diagonal entry 

                receive_buffer_rhs.resize(send_buffer_rhs.size());
                receive_buffer_matrix.resize(send_buffer_matrix.size());

                int other_process = -1;
                if (this_mpi_process == *dof.second.begin())
                    other_process = *std::prev(dof.second.end());
                else
                    other_process = *dof.second.begin();

                MPI_Sendrecv(send_buffer_rhs.data(), send_buffer_rhs.size(), MPI_DOUBLE, other_process, 0,
                            receive_buffer_rhs.data(), receive_buffer_rhs.size(), MPI_DOUBLE, other_process, 0,
                            mpi_communicator, MPI_STATUS_IGNORE);

                MPI_Sendrecv(send_buffer_matrix.data(), send_buffer_matrix.size(), MPI_DOUBLE, other_process, 0,
                            receive_buffer_matrix.data(), receive_buffer_matrix.size(), MPI_DOUBLE, other_process, 0,
                            mpi_communicator, MPI_STATUS_IGNORE);

                for (size_t i = 0; i < receive_buffer_rhs.size(); i += 2)
                {
                    int dof = (int)receive_buffer_rhs[i];
                    double rhs_value = receive_buffer_rhs[i + 1];
                    double matrix_value = receive_buffer_matrix[i];
                    double matrix_value_next = receive_buffer_matrix[i + 1];
                    
                    if (is_red)
                    {
                        system_rhs.at(dof) += rhs_value;
                        system_matrix.add(dof, dof, matrix_value);
                        system_matrix.add(dof, dof + 1 - 2 * (dof <= my_global_dofs[0]), matrix_value_next);
                    }
                    else
                    {
                        system_rhs.erase(dof);
                        system_matrix.erase({dof, dof});
                        system_matrix.erase({dof, dof + 1 - 2 * (dof > my_global_dofs[0])});
                    }
                } 
            }
        }
    }

    const size_t Poisson1D::solve()
    {
        const size_t max_iter = 1000;   // not used atm
        const double tol = 1e-5;
        return solve::parallel::cg(system_matrix, system_rhs, solution, max_iter, tol, mpi_communicator);    
    }

    const double Poisson1D::L2H1norm(const double l2, const double h1)
    {
        
        assert(false);  // not correctly implemented yet

        double val = 0.0;

        std::vector<double> psi_vals(dofs_per_cell);    // container for evaluations of psi
        std::vector<std::vector<double>> dpsi_vals(dofs_per_cell, std::vector<double>(dim, 0.));


        // Loop over all elements
        //for (int k = 0; k < Th.nk; k++) {
        for (auto cell = mesh.cell_begin(); cell != mesh.cell_end(); ++cell) 
        {
            if (cell->get_subdomain() == this_mpi_process) 
            {
                std::vector<size_t> loc2glb(dofs_per_cell);
                for (size_t i = 0; i < dofs_per_cell; ++i) 
                    loc2glb[i] = cell->vertex(i).global_index();

                // Loop over quadrature points
                for (int ipq = 0; ipq < n_quad_pts_per_cell; ++ipq) {

                    const Point<dim> xq(qr[ipq].node);   // quadrature point in reference element

                    if ((int)l2) psi.eval(xq, psi_vals);
                    if ((int)h1) psi.eval_d(*cell, xq, dpsi_vals);

                    const double cint = qr[ipq].weight * cell->get_measure();

                    double uk = 0.;

                    // Loop over dofs
                    for (int i = 0; i < dofs_per_cell; i++) {
                        std::cout << "Process : " << this_mpi_process << " loc2glb[i] : " << loc2glb[i] << " u : " << u(cell->vertex(i)) << " solution : " << solution[loc2glb[i]] << std::endl;    
                        const double diff = solution[loc2glb[i]] - u(cell->vertex(i));
                        const double ul = l2 * psi_vals[i] * diff;
                        uk += ul*ul;

                        for (int j = 0; j < dim; j++) {
                            const double du = h1 * dpsi_vals[i][j] * diff;
                            uk += du*du;
                        }
                    }

                    val += uk * cint;
                }
            }
        }

        // Sum up the contributions from all processes
        double global_val;
        MPI_Allreduce(&val, &global_val, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);

        return std::sqrt(global_val);
    }

    void Poisson1D::output_solution(const std::string &filename) const
    {
        // Let each process write out "solution" to its own file
        std::ofstream file;
        file.open(filename + "_" + std::to_string(this_mpi_process) + ".dat");
        for (size_t i = 0; i < solution.size(); i++)
            file << solution[i] << std::endl;
        file.close();

        std::cout << "Solution written to " << filename + "_" + std::to_string(this_mpi_process) + ".dat" << std::endl;
    }
    
    void Poisson1D::run()
    {
        setup_system();
        assemble_system();
        exchange_shared();      // exchange shared dofs between processes
        const size_t cg_iterations = solve();
        std::cout << "CG iterations: " << cg_iterations << std::endl;
        output_solution("solution");

        // if (this_mpi_process == 0) 
        // {
            
        //     // print soluion
        //     for (size_t i = 0; i < solution.size(); i++)
        //         std::cout << solution[i] << std::endl;
        //     std::cout << std::endl;   
        
        // }

    }

}   // namespace parallel_poisson




