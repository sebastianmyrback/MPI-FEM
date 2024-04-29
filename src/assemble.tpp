#include "quadrature.cpp"

template <typename mesh>
void problem<mesh>::assemble_FEM_matrix(const mesh & Th, const p1_lagrange_1d &psi, const double alpha, const double beta) {

        static const int d = Th.D;
        const int n_quad_pts = 2;

        typedef typename TypeQuadrature<d>::Quadrature Quadrature;
        const Quadrature quadrature;
        typename Quadrature::QuadratureRule qr = quadrature.set_quadrature_rule("trapezoidal");

        const int ndofs = psi.ndof;

        std::vector<double> phi(ndofs);
        std::vector<std::vector<double>> dphi(ndofs, std::vector<double>(0));
        
        // Loop over all elements
        for (int k=0; k<Th.nk; k++) {

            // Get the element
            const typename mesh::elem & K = Th[k];

            // Get the measure of the element
            const double measure = K.measure;

            std::vector<int> loc2glb(ndofs);
            for (int i=0; i<ndofs; i++) {
                loc2glb[i] = K.elem_vertices[i]->idx;
            }

            // Loop over quadrature rule
            for (int ipq = 0; ipq < qr.n; ++ipq) {

                const typename mesh::Rn xq(qr[ipq].node);

                if (alpha) psi.eval(xq, phi);
                if (beta) psi.eval_d(K, xq, dphi);


                for (int i=0; i<ndofs; i++) {
                    for (int j=0; j<ndofs; j++) {
                        
                        

                        mat[std::make_pair(loc2glb[i], loc2glb[j])] += alpha * phi[i] * phi[j] * qr[ipq].weight * measure;
                        
                        // for (int d=0; d<psi.D; d++) {
                        //     mat[std::make_pair(loc2glb[i], loc2glb[j])] += beta * dphi[i][d] * dphi[j][d] * qr[ipq].weight * measure;
                        // }
                    }
                }

            }

        }

    }
