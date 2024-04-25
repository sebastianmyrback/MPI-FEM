#ifndef GC_HPP
#define GC_HPP
#include <assert.h>
#include <iostream>

// Implement the conjugate gradient method
// A is a square matrix of size n x n
// b is a vector of size n
// x is a vector of size n

// A is a symmetric positive definite matrix
// A is stored as a 1D array of size n*n
// A[i*n + j] is the element in the ith row and jth column

// b is stored as a 1D array of size n
// b[i] is the ith element of the vector b

// x is stored as a 1D array of size n
// x[i] is the ith element of the vector x

// The function returns the number of iterations
// The function modifies the vector x

int CG(double *A, double *b, double *x, int n) {
    // Initialize the residual
    double *r = new double[n];
    double *p = new double[n];
    double *Ap = new double[n];
    double alpha, beta, r_dot_r, r_dot_r_new;

    // r = b - A*x
    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        for (int j = 0; j < n; j++) {
            r[i] -= A[i*n + j] * x[j];
        }
        p[i] = r[i];
    }

    r_dot_r = 0.0;
    for (int i = 0; i < n; i++) {
        r_dot_r += r[i] * r[i];
    }

    int max_iter = 1000;
    int iter = 0;
    while (iter < max_iter) {
        // Ap = A*p
        for (int i = 0; i < n; i++) {
            Ap[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Ap[i] += A[i*n + j] * p[j];
            }
        }

        // alpha = r_dot_r / p_dot_Ap
        double p_dot_Ap = 0.0;
        for (int i = 0; i < n; i++) {
            p_dot_Ap += p[i] * Ap[i];
        }
        alpha = r_dot_r / p_dot_Ap;

        // x = x + alpha*p
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
        }

        // r = r - alpha*Ap
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * Ap[i];
        }

        r_dot_r_new = 0.0;
        for (int i = 0; i < n; i++) {
            r_dot_r_new += r[i] * r[i];
        }

        if (r_dot_r_new < 1e-10) {
            break;
        }

        // beta = r_dot_r_new / r_dot_r
        beta = r_dot_r_new / r_dot_r;
        r_dot_r = r_dot_r_new;

        // p = r + beta*p
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        iter++;
    }

    delete[] r;
    delete[] p;
    delete[] Ap;

    return iter;
}


#endif
