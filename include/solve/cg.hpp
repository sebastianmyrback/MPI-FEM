#pragma once

#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include "mpi.h"


namespace solve
{
    // Implement the CG method for solving the linear system Ax = b where A is an std::map<std::pair<int, int>, double>
    // and x, b are std::vector<double>

    namespace serial
    {
        const size_t cg(
            const data_structures::serial::SparseMatrix &A,    // system matrix
            const data_structures::serial::Vector &b,    // right-hand side
            data_structures::serial::Vector &x,          // solution vector
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

            for (int iter = 0; iter < max_iter; iter++) {
                
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

                if (r_dot_r_new < tol) {

                    return iter;
                }

                beta = r_dot_r_new / r_dot_r;

                // Update search direction p
                for (int i = 0; i < n; i++) {
                    p[i] = r[i] + beta * p[i];
                }

                //p = r + beta * p;

            }

            return max_iter;

        }
    }


    namespace parallel
    {
        const size_t cg(
            const data_structures::parallel::SparseMatrix &A_local,    // system matrix
            const data_structures::parallel::Vector &b_local,    // right-hand side
            std::vector<double> &x_local,          // solution vector
            const int max_iter, 
            const double tol) 
        {
        
            // data_structures::parallel::Vector r_local = b_local;      // r0 = b - Ax0 = b
            // data_structures::parallel::Vector p_local = r_local;      // p0 = r0

            std::vector<double> r_local, p_local;
            std::vector<int> loc2glb;

            // loop through b_local which is a map and convert it to a vector
            // assumes the indices are sorted in b_local
            for (auto & [indices, value] : b_local) 
            {
                r_local.push_back(value);
                p_local.push_back(value);
                loc2glb.push_back(indices);
            }

            // // print r_local and its indices
            // for (int i = 0; i < r_local.size(); i++) {
            //     std::cout << r_local[i] << " ";
            // }
            // std::cout << "\n";


            const size_t n_per_process = r_local.size();
            assert(n_per_process == x_local.size());

            int n_total = 0;
            MPI_Allreduce(&n_per_process, &n_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            x_local.assign(n_per_process, 0.0);
            
            std::vector<double> Ap_local;
            double alpha = 0., beta = 0., r_dot_r_local = 0., r_dot_r_new_local = 0., p_dot_Ap_local = 0.;            
            
            for (int iter = 0; iter < max_iter; iter++) 
            {
                
                r_dot_r_local = 0.0;
                for (int i = 0; i < n_per_process; i++) {
                    r_dot_r_local += r_local[i] * r_local[i];
                }

                Ap_local.assign(n_per_process, 0.0);
                for (auto & [indices, value] : A_local) {
                    int row = indices.first;
                    int col = indices.second;

                    // find index in loc2glb that has value col
                    const auto pos_row = std::find(loc2glb.begin(), loc2glb.end(), row);
                    const auto pos_col = std::find(loc2glb.begin(), loc2glb.end(), col);
                    const size_t idx_row_loc = std::distance(loc2glb.begin(), pos_row);
                    const size_t idx_col_loc = std::distance(loc2glb.begin(), pos_col); 
                    

                    // if (it == loc2glb.end()) {
                    //     for (const auto & l : loc2glb) {
                    //         std::cout << l << " ";
                    //     }
                    //     std::cout << "\n";
                    //     std::cout << "it : " << *it << " col : " << col << "\n";
                    //     throw std::invalid_argument("Local and global indices do not match");
                    // }

                    Ap_local[idx_row_loc] += value * p_local[idx_col_loc];

                    //std::cout << "(" << row << ", " << col << ") = " << value << ", idx_row = " << idx_row_loc << ", p_local[" << idx_col_loc << "] = " << p_local[idx_col_loc] << "\n";
                }

                // Gather all elements of p from all processes
                std::vector<double> p_global(n_total);
                MPI_Allgather(p_local.data(), n_per_process, MPI_DOUBLE, p_global.data(), n_per_process, MPI_DOUBLE, MPI_COMM_WORLD);

                // Compute A*p for the local rows of A
                //std::vector<double> Ap_local(n_per_process, 0.0);
                for (int i = 0; i < n_per_process; i++) {
                    for (int j = 0; j < n_total; j++) {
                        Ap_local[i] += A_local[i][j] * p_global[j];
                    }
                }

                p_dot_Ap_local = 0.0;
                for (int i = 0; i < n_per_process; i++) {
                    p_dot_Ap_local += p_local[i] * Ap_local[i];
                }

                double r_dot_r = 0.0, p_dot_Ap = 0.0;
                MPI_Allreduce(&r_dot_r_local, &r_dot_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                alpha = r_dot_r / p_dot_Ap;

                for (int i = 0; i < n_per_process; i++) {
                    x_local[i] += alpha * p_local[i];
                    r_local[i] -= alpha * Ap_local[i];
                }

                r_dot_r_new_local = 0.0;
                for (int i = 0; i < n_per_process; i++) {
                    r_dot_r_new_local += r_local[i] * r_local[i];
                }

                double r_dot_r_new = 0.0;
                MPI_Allreduce(&r_dot_r_new_local, &r_dot_r_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                if (r_dot_r_new < tol) {

                    return iter;
                }

                beta = r_dot_r_new / r_dot_r;

                // Update search direction p
                for (int i = 0; i < n_per_process; i++) {
                    p_local[i] = r_local[i] + beta * p_local[i];
                }


            }

            return max_iter;

        }
    }


}
