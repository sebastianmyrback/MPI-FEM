#include <iostream>
#include "problem.hpp"
#include "export.hpp"
#include "cg.hpp"
#include "norm.hpp"
//#include "matplotlibcpp.h"
//namespace plt = matplotlibcpp;

// define a function to be used as a source term f(x) = 8*pi^2*sin(2*pi*x)
double f(const Rd<1> & x) {
    //return 8.0 * M_PI * M_PI * sin(2.0 * M_PI * x[0]);
    return (M_PI*M_PI*(361*cos((19*M_PI*x[0])/10) - 441*cos((21*M_PI*x[0])/10)))/100;
}

double u(const Rd<1> & x) {
    //return 2*sin(2.0 * M_PI * x[0]);
    return 2*sin(2*M_PI*x[0])*sin(M_PI*x[0]/10);
}

int main() {
    
    // Create a mesh object
    const int n_refinements = 7;
    int n = 7;     // number of elements
    const int n_threads = 1;
    const double a = 0., b = 5.;

    std::vector<double> l2_errors(n_refinements, 0.), h1_errors(n_refinements, 0.);

    std::vector<double> mesh_vertices, mesh_sizes, mesh_sizes2, uh, uexact, diff;

    QuadratureRule<1> midpoint(1, {QuadraturePoint<1>(Rd<1>({0.5}), 1.)});
    
    QuadratureRule<1> trapezoidal(2, {QuadraturePoint<1>(Rd<1>({0.0}), 0.5), QuadraturePoint<1>(Rd<1>({1.0}), 0.5)});

    QuadratureRule<1> simpson(3, {QuadraturePoint<1>(Rd<1>({0.0}), 1./3), QuadraturePoint<1>(Rd<1>({.5}), 1./3), QuadraturePoint<1>(Rd<1>({1.}), 1./3)});
    
    for (int i = 0; i < n_refinements; i++) {

        std::cout << i + 1 << " / " << n_refinements << std::endl;
        
        lagrange_1d<1> psi;
        
        mesh1d Th(a, b, n);

        problem<mesh1d> prob(n_threads, std::make_shared<mesh1d>(Th));

        mesh_sizes.push_back(10*Th.h);
        mesh_sizes2.push_back(10*Th.h*Th.h);

        std::vector<double> mesh_vertices, uh, uexact, diff;
        for (int i = 0; i < Th.nv; i++) {
            mesh_vertices.push_back(Th(i).x[0]);
            uexact.push_back(u(Th(i).x));
        }

        std::vector<int> dirichlet_lbs = {1, 2};

        prob.assemble_rhs(midpoint, psi, f);

        //! SOMETHING IS PROBABLY WRONG WITH THE QUADRATURE RULES, SINCE MIDPOINT IS TWICE AS GOOD 
        //! AS SIMPSON AND TRAPEZOIDAL

        // assemble stiffness matrix
        prob.assemble_FEM_matrix(midpoint, psi, 0., 1., dirichlet_lbs);
        //prob.assemble_FEM_matrix(Th, trapezoidal, psi, 0., 1.);
        
        //prob.set_dirichlet(Th, u);      // currently only works for homogeneous Dirichlet boundary conditions

        const double tol = 1e-10;
        const int max_iter = 1000;

        // Chose initial guess u0 as zeros and the exact boundary values
        std::vector<double> u0(Th.nv, 0.);  
        // for (int i = 0; i < Th.nv; i++) {
        //     const vertex<1> & v(Th(i));
        //     int lab = v.vertex_label;
        //     if (lab != 0) {
        //         u0[i] = u(v.x);
        //     }
        // }

        // Solve the linear system using the conjugate gradient method
        uh = cg(prob.mat, prob.rhs, u0, max_iter, tol);

        // Compute the difference between the exact and approximate solution coefficients
        diff.clear();
        for (int i = 0; i < Th.nv; i++) {
            diff.push_back(uexact[i] - uh[i]);
        }

        // Compute the L2 and H1 errors
        l2_errors[i] = L2H1norm(Th, trapezoidal, psi, diff, 1., 0.);
        h1_errors[i] = L2H1norm(Th, trapezoidal, psi, diff, 0., 1.);

        n *= 2;

        matlab::save(prob.mat, "matrix.dat");
        matlab::save(prob.rhs, "rhs.dat");
        matlab::save(mesh_vertices, "x.dat");
        matlab::save(uh, "uh.dat");

    }

    // print the errors
    std::cout << "L2 errors: ";
    for (auto & e : l2_errors) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

    std::cout << "H1 errors: ";
    for (auto & e : h1_errors) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

    // std::cout << "Mesh vertices: ";
    // for (auto & x : mesh_vertices) {
    //     std::cout << x << " ";
    // }  
    // std::cout << std::endl;



    // plt::plot(mesh_vertices, uh, "*", {{"label", "uh"}});
    // plt::plot(mesh_vertices, uexact, {{"label", "u exact"}});
    // plt::legend();
    // plt::show();

    // plt::loglog(mesh_sizes, mesh_sizes2, {{"label", "h^2"}});
    // plt::loglog(mesh_sizes, mesh_sizes, {{"label", "h"}});
    // plt::loglog(mesh_sizes, l2_errors, "*", {{"label", "L2 error"}});
    // plt::loglog(mesh_sizes, h1_errors, "^", {{"label", "H1 error"}});
    // plt::legend();
    // plt::show();

    return 0;
}