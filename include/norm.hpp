#ifndef NORM_HPP
#define NORM_HPP

#include "basis_functions.hpp"
#include "quadrature.hpp"

template <typename mesh, int d, int degree>
static double L2H1norm(const mesh &Th, const QuadratureRule<d> &qr, const basis_function<d, degree> & psi, const std::vector<double> &u, const double l2, const double h1) {

    typedef typename mesh::Rn Rn;

    double val = 0.0;

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int ndofs = psi.ndof;             // number of dofs per element
    std::vector<double> psi_vals(ndofs);    // container for evaluations of psi
    std::vector<std::vector<double>> dpsi_vals(ndofs, std::vector<double>(dim, 0.));    // ndofs x dim matrix


    // Loop over all elements
    //for (int k = 0; k < Th.nk; k++) {
    for (auto cell = Th.cell_begin(); cell != Th.cell_end(); ++cell) 
    {

        // Create map from local to global dofs
        std::vector<int> loc2glb(ndofs);
        for (int i=0; i<ndofs; i++)
            loc2glb[i] = cell->vertex(i).global_index();

        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

            const Rn xq(qr[ipq].node);   // quadrature point in reference element

            if (l2) psi.eval(xq, psi_vals);
            if (h1) psi.eval_d(*cell, xq, dpsi_vals);

            const double cint = qr[ipq].weight * cell->get_measure();

            double uk = 0.;

            // Loop over dofs
            for (int i = 0; i < ndofs; i++) {
                const double ul = l2 * psi_vals[i] * u[loc2glb[i]];
                uk += ul*ul;

                for (int j = 0; j < dim; j++) {
                    const double du = h1 * dpsi_vals[i][j] * u[loc2glb[i]];
                    uk += du*du;
                }
            }

            val += uk * cint;
        }
    }

    return std::sqrt(val);

}

#endif // NORM_HPP