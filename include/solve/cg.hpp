#pragma once

#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include "mpi.h"
#include "../utilities/mpi_util.hpp"


namespace solve
{
    // Implement the CG method for solving the linear system Ax = b where A is an std::map<std::pair<int, int>, double>
    // and x, b are std::vector<double>

    namespace serial
    {
        const size_t cg(
            const data_structures::serial::SparseMatrix &A,    // system matrix
            const data_structures::serial::Vector &b,    // right-hand side
            std::vector<double> &x,          // solution vector
            const int max_iter, 
            const double tol) 
        {
        
            const int n = b.size();

            if (n != x.size()) {
                throw std::invalid_argument("Mismatch in sizes of x and b");
            }

            //std::vector<double> x(u0);  // x0 = 0
            data_structures::serial::Vector r = b;      // r0 = b - Ax0 = b
            data_structures::serial::Vector p = r;      // p0 = r0
            data_structures::serial::Vector Ap(n, 0.0);
            double alpha = 0., beta = 0., r_dot_r = 0., r_dot_r_new = 0., p_dot_Ap = 0.;

            int iter = 0;
            for (iter = 0; iter < n; iter++) {
                
                // Compute r^T*r
                r_dot_r = 0.0;
                for (int i = 0; i < n; i++) {
                    r_dot_r += r[i] * r[i];
                }

                //r_dot_r = r.dot(r);

                // Compute A*pk
                Ap.assign(n, 0.0);
                for (auto & [indices, value] : A) {
                    int row = indices.first;
                    int col = indices.second;

                    Ap[row] += value * p[col];
                }
                //Ap = A * p;

                // Compute pk^T*(A*pk)
                p_dot_Ap = 0.0;
                for (int i = 0; i < n; i++) {
                    p_dot_Ap += p[i] * Ap[i];
                }

                //p_dot_Ap = p.dot(Ap);

                alpha = r_dot_r / p_dot_Ap;


                // Update x and r
                for (int i = 0; i < n; i++) {
                    x[i] += alpha * p[i];
                    r[i] -= alpha * Ap[i];
                }
                // x += alpha * p;
                // r -= alpha * Ap;

                r_dot_r_new = 0.0;
                for (int i = 0; i < n; i++) {
                    r_dot_r_new += r[i] * r[i];
                }
                
                //r_dot_r_new = r.dot(r);

                if (std::sqrt(r_dot_r_new) < tol) {
		            std::cout << "Residual after " << iter << " iterations: " << std::sqrt(r_dot_r_new) << "\n";

                    return iter;
                }

                beta = r_dot_r_new / r_dot_r;

                // Update search direction p
                for (int i = 0; i < n; i++) {
                    p[i] = r[i] + beta * p[i];
                }

                //p = r + beta * p;

            }

	        std::cout << "Residual after " << iter-1 << " iterations (max): " << std::sqrt(r_dot_r_new) << "\n";

            return iter;

        }
    }


   namespace parallel
    {
        const size_t cg(
            const data_structures::parallel::SparseMatrix &A_local,    // system matrix
            const data_structures::parallel::Vector &b_local,          // right-hand side
            std::vector<double> &x_local,                              // solution vector 
            const int max_iter, 
            const double tol,
            MPI_Comm mpi_communicator) 
        {
            const int this_rank = mpi_util::this_mpi_process(mpi_communicator);
            const int n_processes = mpi_util::n_mpi_processes(mpi_communicator);

            std::vector<double> r_local, p_local;
            r_local.reserve(b_local.size());
            p_local.reserve(b_local.size());

            std::unordered_map<int, size_t> loc2glb;
            size_t i = 0;
            for (auto & [indices, value] : b_local) 
            {
                r_local.push_back(value);
                p_local.push_back(value);
                loc2glb[indices] = i++;
            }

            const size_t n_local = r_local.size();  // number of rows for current process
            x_local.assign(n_local, 0.);            // initial guess always zero

            int n_total = 0;
            MPI_Allreduce(&n_local, &n_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            // Gather the number of elements per process and their relative position to the global vector
            std::vector<int> recvcounts(n_processes), displs(n_processes);
            MPI_Allgather(&n_local, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
            displs[0] = 0;
            for (int i = 1; i < n_processes; ++i) {
                displs[i] = displs[i - 1] + recvcounts[i - 1];
            }
            
            std::vector<double> Ap_local(n_local), p_global(n_total);
            double alpha = 0., beta = 0., r_dot_r_local = 0., r_dot_r_new_local = 0., p_dot_Ap_local = 0., r_dot_r_new = 0.;            

            int iter = 0;
            for (iter = 0; iter < n_total; iter++)
            {
                r_dot_r_local = 0.0;
                for (int i = 0; i < n_local; i++) {
                    r_dot_r_local += r_local[i] * r_local[i];
                }

                MPI_Allgatherv(p_local.data(), n_local, MPI_DOUBLE, p_global.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

                Ap_local.assign(n_local, 0.0);

                for (auto & [indices, value] : A_local) 
                {
                    int row = indices.first;
                    int col = indices.second;

                    const auto pos_row_it = loc2glb.find(row);
                    if (pos_row_it != loc2glb.end()) {
                        const size_t idx_row_loc = pos_row_it->second;
                        Ap_local[idx_row_loc] += value * p_global[col];
                    }
                    
                }

                p_dot_Ap_local = 0.0;
                for (int i = 0; i < n_local; i++) {
                    p_dot_Ap_local += p_local[i] * Ap_local[i];
                }

                double r_dot_r = 0.0, p_dot_Ap = 0.0;
                double local_values[2] = {r_dot_r_local, p_dot_Ap_local};
                double global_values[2] = {0.0, 0.0};
                MPI_Allreduce(local_values, global_values, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                r_dot_r = global_values[0];
                p_dot_Ap = global_values[1];

                alpha = r_dot_r / p_dot_Ap;

                for (int i = 0; i < n_local; i++) {
                    x_local[i] += alpha * p_local[i];
                    r_local[i] -= alpha * Ap_local[i];
                }

                r_dot_r_new_local = 0.0;
                for (int i = 0; i < n_local; i++) {
                    r_dot_r_new_local += r_local[i] * r_local[i];
                }

                r_dot_r_new = 0.0;
                MPI_Allreduce(&r_dot_r_new_local, &r_dot_r_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // if (this_rank == 0)
                //     std::cout << "Iteration " << iter << ": Residual: " << std::sqrt(r_dot_r_new) << "\n";

                if (std::sqrt(r_dot_r_new) < tol) {
		            if (this_rank == 0)
            		    std::cout << "Final residual after " << iter << " iterations: " << std::sqrt(r_dot_r_new) << "\n";
		            return iter;
                }

                beta = r_dot_r_new / r_dot_r;

                for (int i = 0; i < n_local; i++) {
                    p_local[i] = r_local[i] + beta * p_local[i];
                }

            }

	    if (this_rank == 0)
	        std::cout << "Final residual after " << iter-1 << " iterations (max): " << std::sqrt(r_dot_r_new) << "\n";
            return iter;
        }
    }

}
