#include "quadrature.cpp"

template <typename mesh>
template <int d, int deg>
void problem<mesh>::assemble_FEM_matrix(const mesh & Th, const QuadratureRule<d> &qr, const basis_function<d, deg> &psi, const double mass, const double stiffness) {

    // clear the matrix
    mat.clear();

    // alpha = 1 -> mass matrix
    // beta = 1 -> stiffness matrix

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int ndofs = psi.ndof;

    std::vector<double> psi_vals(ndofs);
    std::vector<std::vector<double>> dpsi_vals(ndofs, std::vector<double>(dim, 0.));    // ndofs x dim matrix
    
    // Loop over all elements
    for (int k = 0; k < Th.nk; k++) {

        // Get the element
        const elem & K = Th[k];

        // Get the measure of the element
        const double measure = K.measure;

        // Create map from local to global dofs
        std::vector<int> loc2glb(ndofs);
        for (int i=0; i<ndofs; i++)
            loc2glb[i] = K.elem_vertices[i]->glb_idx;
        
        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

            const Rn xq(qr[ipq].node);   // quadrature point in reference element

            if (mass) psi.eval(xq, psi_vals);
            if (stiffness) psi.eval_d(K, xq, dpsi_vals);

            // Loop over dofs
            for (int i = 0; i < ndofs; i++) {
                for (int j = 0; j < ndofs; j++) {
                    
                    if (mass)
                        mat[std::make_pair(loc2glb[i], loc2glb[j])] += mass * psi_vals[i] * psi_vals[j] * qr[ipq].weight * measure;
                    
                    if (stiffness) {
                        for (int dm = 0; dm < dim; dm++) {
                            mat[std::make_pair(loc2glb[i], loc2glb[j])] += stiffness * dpsi_vals[i][dm] * dpsi_vals[j][dm] * qr[ipq].weight * measure;
                        }
                    }
                    
                }
            }

        }

    }

}



template <typename mesh>
template <int d, int deg>
void problem<mesh>::assemble_rhs(const mesh & Th, const QuadratureRule<d> &qr, const basis_function<d, deg> &psi, const std::function<double(const typename mesh::Rn &)> & f) {

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int ndofs = psi.ndof;             // number of dofs per element

    std::vector<double> psi_vals(ndofs);    // container for evaluations of psi

    rhs.resize(Th.nv, 0.);                 // resize rhs to number of dofs

    // Loop over all elements
    for (int k = 0; k < Th.nk; k++) {

        // Get the current element
        const elem & K = Th[k];

        // Get the measure of the element
        const double measure = K.measure;

        // Create map from local to global dofs
        std::vector<int> loc2glb(ndofs);
        for (int i=0; i<ndofs; i++)
            loc2glb[i] = K.elem_vertices[i]->glb_idx;
        
        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

            const Rn xq(qr[ipq].node);   // quadrature point in reference element

            Rn x;
            K(xq, x);           // map xq to quadrature point in physical element x

            psi.eval(xq, psi_vals);                     // evaluate psi at xq -> store in psi_vals
            
            // Loop over dofs
            for (int i = 0; i < ndofs; i++) {
                // std::cout << "x = " << x << std::endl;
                // std::cout << "f(x) = " << f(x) << std::endl;
                // getchar();
                rhs[loc2glb[i]] += psi_vals[i] * f(x) * qr[ipq].weight * measure;
            }

        }

    }

}



template <typename mesh>
void problem<mesh>::set_dirichlet(const mesh & Th, const std::function<double(const typename mesh::Rn &)> & g) {

    for (int i = 0; i < Th.nv; i++) {
        if (Th.mesh_vertices[i].vertex_label != 0) {
            rhs[i] = g(Th.mesh_vertices[i].x);

            std::cout << rhs[i] << std::endl;

            // Zero out the row and column corresponding to the Dirichlet boundary
            for (auto & [key, value] : mat) {
                if (key.first == i || key.second == i) {
                    value = 0.0;
                }
            }
            mat[std::make_pair(i, i)] = 1.0;
        }
    }

}


