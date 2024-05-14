#pragma once

#include "../utilities/data_structures.hpp"

namespace utilities 
{
    template <int d>
    struct DirichletBC {
        const double (*g)(const Point<d> &);    // Dirichlet boundary function
        std::vector<int> lbs;                   // labels of Dirichlet boundaries
        bool set_dirichlet = false;             // if false, no Dirichlet boundary condition is set strongly
    };

    void distribute_local_to_global(
        data_structures::serial::DenseMatrix &Ak,
        data_structures::serial::Vector &fk,
        data_structures::serial::SparseMatrix &system_matrix,
        data_structures::serial::Vector &system_rhs,
        const std::vector<size_t> &loc2glb,
        const std::map<int, double> &boundary_data)
    {

        if (loc2glb.size() != Ak.size() || loc2glb.size() != fk.size())
            throw std::invalid_argument("Mismatch in sizes of loc2glb, Ak, and fk");


        for (int i = 0; i < loc2glb.size(); i++)
        {

            // Add local rhs vector to global rhs vector
            system_rhs(loc2glb[i]) += fk[i];

            
            // Check if dof i is a Dirichlet dof
            double avg_diag = 0.; 
            bool is_row_dof_dirichlet = false;
            is_row_dof_dirichlet = (boundary_data.find(loc2glb[i]) != boundary_data.end()); 
        

            for (int j = 0; j < loc2glb.size(); j++)
            {
            
                avg_diag += Ak(j, j) / loc2glb.size();

                // Check if dof j is a Dirichlet dof
                bool is_col_dof_dirichlet = false;
                is_col_dof_dirichlet = (boundary_data.find(loc2glb[j]) != boundary_data.end());

                // If i or j is a Dirichlet dof:
                // set Dirichlet condition in local matrix and
                // change the corresponding entry in the rhs vector
                if (is_row_dof_dirichlet || is_col_dof_dirichlet) 
                {
                    // Subtract dirichlet column of A times boundary value from rhs
                    if (is_col_dof_dirichlet) 
                    {
                        system_rhs(loc2glb[i]) -= Ak(i, j) * boundary_data.at(loc2glb[j]);
                    }

                    // Set local row and column to zero
                    Ak(i, j) = 0.0;

                    // Set diagonal entry to average of all other local diagonal entries
                    // modify rhs accordingly
                    if (i == j) 
                    {
                        Ak(i, j) = avg_diag;
                        system_rhs(loc2glb[i]) = avg_diag * boundary_data.at(loc2glb[j]);
                    }
                }

                // Add local matrix entry to global matrix
                //system_matrix[std::make_pair(loc2glb[i], loc2glb[j])] += Ak(i, j);
                system_matrix.add(loc2glb[i], loc2glb[j], Ak(i, j));
                
            }
        }
    }


    namespace parallel_dofs 
    {
        void distribute_local_to_global(
            data_structures::serial::DenseMatrix &Ak,
            data_structures::serial::Vector &fk,
            data_structures::parallel::SparseMatrix &system_matrix,
            data_structures::parallel::Vector &system_rhs,
            const std::vector<size_t> &loc2glb,
            const std::map<int, double> &boundary_data
            //,
            //MPI_Comm mpi_communicator
            )
        {

            if (loc2glb.size() != Ak.size() || loc2glb.size() != fk.size())
                throw std::invalid_argument("Mismatch in sizes of loc2glb, Ak, and fk");


            for (int i = 0; i < loc2glb.size(); i++)
            {

                // Add local rhs vector to global rhs vector
                system_rhs(loc2glb[i]) += fk[i];

                
                // Check if dof i is a Dirichlet dof
                double avg_diag = 0.; 
                bool is_row_dof_dirichlet = false;
                is_row_dof_dirichlet = (boundary_data.find(loc2glb[i]) != boundary_data.end()); 
            

                for (int j = 0; j < loc2glb.size(); j++)
                {
                
                    avg_diag += Ak(j, j) / loc2glb.size();

                    // Check if dof j is a Dirichlet dof
                    bool is_col_dof_dirichlet = false;
                    is_col_dof_dirichlet = (boundary_data.find(loc2glb[j]) != boundary_data.end());

                    // If i or j is a Dirichlet dof:
                    // set Dirichlet condition in local matrix and
                    // change the corresponding entry in the rhs vector
                    if (is_row_dof_dirichlet || is_col_dof_dirichlet) 
                    {
                        // Subtract dirichlet column of A times boundary value from rhs
                        if (is_col_dof_dirichlet) 
                        {
                            system_rhs(loc2glb[i]) -= Ak(i, j) * boundary_data.at(loc2glb[j]);
                        }

                        // Set local row and column to zero
                        Ak(i, j) = 0.0;

                        // Set diagonal entry to average of all other local diagonal entries
                        // modify rhs accordingly
                        if (i == j) 
                        {
                            Ak(i, j) = avg_diag;
                            system_rhs(loc2glb[i]) = avg_diag * boundary_data.at(loc2glb[j]);
                        }
                    }

                    // Add local matrix entry to global matrix
                    //system_matrix[std::make_pair(loc2glb[i], loc2glb[j])] += Ak(i, j);
                    system_matrix.add(loc2glb[i], loc2glb[j], Ak(i, j));
                    
                }
            }
        }
    }


}