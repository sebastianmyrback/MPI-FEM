#include <iostream>
#include "fe/fem.hpp"
#include "utilities/export.hpp"
#include "utilities/norm.hpp"
#include "solve/cg.hpp"
#include <chrono>

// #include "matplotlibcpp.h"
// namespace plt = matplotlibcpp;




// define a function to be used as a source term f(x) = 8*pi^2*sin(2*pi*x)
const double f(const Point<1> & x) {
    //return 8.0 * M_PI * M_PI * sin(2.0 * M_PI * x[0]);
    return (M_PI*M_PI*(361*cos((19*M_PI*x[0])/10) - 441*cos((21*M_PI*x[0])/10)))/100;
}

const double u(const Point<1> & x) {
    //return 2*sin(2.0 * M_PI * x[0]);
    return 2*sin(2*M_PI*x[0])*sin(M_PI*x[0]/10) + 10;
}


int main() {

    using namespace quadrature;
    auto start = std::chrono::high_resolution_clock::now();

    // Create a mesh object
    const int n_refinements = 7;
    int n = 10;     // number of elements
    const int n_threads = 1;
    const double a = -0.63, b = 5.27;
    //const double a = 0., b = 1.;
    
    std::vector<double> l2_errors(n_refinements, 0.), h1_errors(n_refinements, 0.), mesh_sizes(n_refinements, 0.), mesh_sizes_sq(n_refinements, 0.);
    data_structures::serial::Vector mesh_vertices, uh, uexact, diff;

    for (int i = 0; i < n_refinements; i++) {

        std::cout << i + 1 << " / " << n_refinements << std::endl;
        
        P1Lagrange1D<1> psi;
        const Mesh1D Th(a, b, n);
        //FEM<Mesh1D> prob(n_threads, std::make_shared<Mesh1D>(Th));
        FEM<Mesh1D> prob(n_threads, &Th);

        mesh_sizes[i] = Th.get_h();
        
        mesh_vertices.assign(Th.get_nverts(), 0.);
        uh.assign(Th.get_nverts(), 0.);
        uexact.assign(Th.get_nverts(), 0.);
        diff.assign(Th.get_nverts(), 0.);
        

        for (int v = 0; v < Th.get_nverts(); v++) {
            mesh_vertices[v] = Th.get_vertices()[v][0];
            uexact[v] = u(Th.get_vertices()[v]);
        }

        utilities::DirichletBC<1> bc;
        bc.g = u;
        bc.lbs = {1, 2};
        bc.set_dirichlet = true;

        prob.assemble_stiffness_system(midpoint, psi, f, bc);

        const double tol = 1e-10;
        const int max_iter = 1000;

        const int cg_iterations = solve::serial::cg(prob.mat, prob.rhs, uh, max_iter, tol);

        diff = uexact - uh;

        std::cout << "CG iterations: " << cg_iterations << std::endl;

        // Compute the L2 and H1 errors
        l2_errors[i] = L2H1norm(Th, gauss_lobatto6, psi, diff, 1., 0.);
        h1_errors[i] = L2H1norm(Th, gauss_lobatto6, psi, diff, 0., 1.);

        n *= 2;

        //gnuplot::write_cells(Th, "Th.dat");

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


    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Elapsed time: " <<  duration.count() << " [ms]" << std::endl;


    // plt::plot(mesh_vertices, uh, "*", {{"label", "uh"}});
    // plt::plot(mesh_vertices, uexact, {{"label", "u exact"}});
    // plt::legend();
    // plt::show();

    // plt::loglog(mesh_sizes, mesh_sizes_sq, {{"label", "h^2"}});
    // plt::loglog(mesh_sizes, mesh_sizes, {{"label", "h"}});
    // plt::loglog(mesh_sizes, l2_errors, "*", {{"label", "L2 error"}});
    // plt::loglog(mesh_sizes, h1_errors, "^", {{"label", "H1 error"}});
    // plt::legend();
    // plt::show();


    return 0;
}