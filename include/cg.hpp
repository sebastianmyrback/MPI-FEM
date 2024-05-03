#ifndef GC_HPP
#define GC_HPP
#include <assert.h>
#include <iostream>
#include <vector>
#include <map>

// Implement the CG method for solving the linear system Ax = b where A is an std::map<std::pair<int, int>, double>
// and x, b are std::vector<double>
std::vector<double> cg(const std::map<std::pair<int, int>, double> &A, const std::vector<double> &b, const int max_iter, const double tol) {
    
    const int n = b.size();

    std::vector<double> x(n, 0.0);  // x0 = 0
    std::vector<double> r = b;      // r0 = b - Ax0 = b
    std::vector<double> p = r;      // p0 = r0
    std::vector<double> Ap(n, 0.0);
    double alpha = 0., beta = 0., r_dot_r = 0., r_dot_r_new = 0., p_dot_Ap = 0.;

    for (int k = 0; k < max_iter; k++) {
        
        // Compute r^T*r
        r_dot_r = 0.0;
        for (int i = 0; i < n; i++) {
            r_dot_r += r[i] * r[i];
        }

        // Compute A*pk
        Ap.assign(n, 0.0);
        for (auto & [indices, value] : A) {
            int row = indices.first;
            int col = indices.second;

            Ap[row] += value * p[col];
        }

        // Compute pk^T*(A*pk)
        p_dot_Ap = 0.0;
        for (int i = 0; i < n; i++) {
            p_dot_Ap += p[i] * Ap[i];
        }

        alpha = r_dot_r / p_dot_Ap;


        // Update x and r
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        r_dot_r_new = 0.0;
        for (int i = 0; i < n; i++) {
            r_dot_r_new += r[i] * r[i];
        }

        if (r_dot_r_new < tol) {
            std::cout << "CG converged after " << k + 1 << " iterations\n";
            return x;
        }

        beta = r_dot_r_new / r_dot_r;

        // Update search direction p
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

    }

    return x;

}
#endif
