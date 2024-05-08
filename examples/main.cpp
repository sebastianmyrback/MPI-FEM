#include <iostream>
#include "problem.hpp"
#include "export.hpp"
#include "cg.hpp"
#include "norm.hpp"
#include <chrono>


#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

const QuadratureRule<1> midpoint(1, {
    QuadraturePoint<1>(Rd<1>({0.5}), 1.)
    });

const QuadratureRule<1> trapezoidal(2, {
    QuadraturePoint<1>(Rd<1>({0.0}), 0.5), 
    QuadraturePoint<1>(Rd<1>({1.0}), 0.5)
    });

const QuadratureRule<1> simpson(3, {
    QuadraturePoint<1>(Rd<1>({0.0}), 1./3), 
    QuadraturePoint<1>(Rd<1>({0.5}), 1./3), 
    QuadraturePoint<1>(Rd<1>({1.0}), 1./3)
    });

const QuadratureRule<1> gauss_lobatto6(6, {
    QuadraturePoint<1>(Rd<1>({0.0}), 0.03333333333333333), 
    QuadraturePoint<1>(Rd<1>({0.11747233803526763}), 0.1892374781489235), 
    QuadraturePoint<1>(Rd<1>({0.3573842417596774}), 0.2774291885177432),
    QuadraturePoint<1>(Rd<1>({0.6426157582403226}), 0.2774291885177432), 
    QuadraturePoint<1>(Rd<1>({0.8825276619647324}), 0.1892374781489235), 
    QuadraturePoint<1>(Rd<1>({1.0}), 0.03333333333333333)
    });




// define a function to be used as a source term f(x) = 8*pi^2*sin(2*pi*x)
double f(const Rd<1> & x) {
    //return 8.0 * M_PI * M_PI * sin(2.0 * M_PI * x[0]);
    return (M_PI*M_PI*(361*cos((19*M_PI*x[0])/10) - 441*cos((21*M_PI*x[0])/10)))/100;
}

double u(const Rd<1> & x) {
    //return 2*sin(2.0 * M_PI * x[0]);
    return 2*sin(2*M_PI*x[0])*sin(M_PI*x[0]/10) + 10;
}


int main() {

    auto start = std::chrono::high_resolution_clock::now();

    // Create a mesh object
    const int n_refinements = 7;
    int n = 10;     // number of elements
    const int n_threads = 1;
    const double a = -0.63, b = 5.27;
    

    std::vector<double> l2_errors(n_refinements, 0.), h1_errors(n_refinements, 0.), mesh_sizes(n_refinements, 0.), mesh_sizes_sq(n_refinements, 0.);
    std::vector<double> mesh_vertices, uh, uexact, diff;


    for (int i = 0; i < n_refinements; i++) {

        std::cout << i + 1 << " / " << n_refinements << std::endl;
        
        lagrange_1d<1> psi;
        
        Mesh1D Th(a, b, n);

        problem<Mesh1D> prob(n_threads, std::make_shared<Mesh1D>(Th));

        mesh_sizes[i] = Th.get_h();

        mesh_vertices.clear();
        uh.clear();
        uexact.clear();
        diff.clear();

        for (auto vertex = Th.vertex_begin(); vertex != Th.vertex_end(); ++vertex) {
            mesh_vertices.push_back((*vertex)[0]);
            uexact.push_back(u((*vertex)[0]));
        }

        dirichlet_bc<1> bc;
        bc.g = u;
        bc.lbs = {1, 2};
        bc.set_dirichlet = true;

        const double mass = 0;
        const double stiffness = 1;

        prob.assemble_rhs(midpoint, psi, f);
        prob.assemble_FEM_matrix(midpoint, psi, mass, stiffness, bc);
        //prob.assemble_FEM_matrix(trapezoidal, psi, mass, stiffness);
        
        const double tol = 1e-10;
        const int max_iter = 1000;

        // Chose initial guess u0 as zeros
        std::vector<double> u0(Th.get_nverts(), 0.);  
        
        // Solve the linear system using the conjugate gradient method
        uh = cg(prob.mat, prob.rhs, u0, max_iter, tol);

        // Compute the difference between the exact and approximate solution coefficients
        diff.clear();
        for (int i = 0; i < Th.get_nverts(); i++) {
            diff.push_back(uexact[i] - uh[i]);
        }

        // Compute the L2 and H1 errors
        l2_errors[i] = L2H1norm(Th, gauss_lobatto6, psi, diff, 1., 0.);
        h1_errors[i] = L2H1norm(Th, gauss_lobatto6, psi, diff, 0., 1.);

        n *= 2;

        matlab::save(prob.mat, "matrix.dat");
        matlab::save(prob.rhs, "rhs.dat");
        matlab::save(mesh_vertices, "x.dat");
        matlab::save(uh, "uh.dat");

    }

    std::cout << "Mesh sizes: ";
    for (int i = 0; i < n_refinements; i++) {
        std::cout << mesh_sizes[i] << " ";

        mesh_sizes_sq[i] = 10 * mesh_sizes[i] * mesh_sizes[i];
        mesh_sizes[i] *= 10;
    }
    std::cout << std::endl;

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



    plt::plot(mesh_vertices, uh, "*", {{"label", "uh"}});
    plt::plot(mesh_vertices, uexact, {{"label", "u exact"}});
    plt::legend();
    plt::show();

    plt::loglog(mesh_sizes, mesh_sizes_sq, {{"label", "h^2"}});
    plt::loglog(mesh_sizes, mesh_sizes, {{"label", "h"}});
    plt::loglog(mesh_sizes, l2_errors, "*", {{"label", "L2 error"}});
    plt::loglog(mesh_sizes, h1_errors, "^", {{"label", "H1 error"}});
    plt::legend();
    plt::show();

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Elapsed time: " <<  duration.count() << " [ms]" << std::endl;

    return 0;
}