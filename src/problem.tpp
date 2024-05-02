#include "quadrature.cpp"

template <typename mesh>
template <int d, int deg>
void problem<mesh>::assemble_FEM_matrix(const mesh & Th, const QuadratureRule<d> &qr, const basis_function<d, deg> &psi, const double alpha, const double beta) {

    // alpha = 1 -> mass matrix
    // beta = 1 -> stiffness matrix

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int ndofs = psi.ndof;

    std::vector<double> phi(ndofs);
    std::vector<std::vector<double>> dphi(ndofs, std::vector<double>(dim, 0.));    // ndofs x dim matrix
    
    // Loop over all elements
    for (int k = 0; k < Th.nk; k++) {

        // Get the element
        const typename mesh::elem & K = Th[k];

        // Get the measure of the element
        const double measure = K.measure;

        // Create map from local to global dofs
        std::vector<int> loc2glb(ndofs);
        for (int i=0; i<ndofs; i++)
            loc2glb[i] = K.elem_vertices[i]->glb_idx;
        
        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

            const typename mesh::Rn xq(qr[ipq].node);   // quadrature point in reference element

            if (alpha) psi.eval(xq, phi);
            if (beta) psi.eval_d(K, xq, dphi);

            // Loop over dofs
            for (int i = 0; i < ndofs; i++) {
                for (int j = 0; j < ndofs; j++) {
                    
                    if (alpha)
                        mat[std::make_pair(loc2glb[i], loc2glb[j])] += alpha * phi[i] * phi[j] * qr[ipq].weight * measure;
                    
                    if (beta) {
                        for (int dm = 0; dm < dim; dm++) {
                            mat[std::make_pair(loc2glb[i], loc2glb[j])] += beta * dphi[i][dm] * dphi[j][dm] * qr[ipq].weight * measure;
                        }
                    }
                    
                }
            }

        }

    }

}
