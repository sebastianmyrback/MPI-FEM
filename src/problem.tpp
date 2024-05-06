#include "quadrature.cpp"

template <typename mesh>
template <int d, int deg>
void problem<mesh>::assemble_FEM_matrix(const QuadratureRule<d> &qr, const basis_function<d, deg> &psi, const double mass, const double stiffness, const std::vector<int> &dirichlet_lbls) {

    // clear the matrix
    mat.clear();

    // alpha = 1 -> mass matrix
    // beta = 1 -> stiffness matrix

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int dofs_per_elem = psi.ndof;

    std::vector<double> psi_vals(dofs_per_elem);
    std::vector<std::vector<double>> dpsi_vals(dofs_per_elem, std::vector<double>(dim, 0.));    // ndofs x dim matrix
    
    // Loop over all elements
    for (int k = 0; k < Th->nk; k++) {

        double avg_diag = 0.;    // average of diagonal entries

        // Get the element
        const elem & K = (*Th)[k];

        // Get the measure of the element
        const double measure = K.measure;

        // Create map from local to global dofs and corresponding vertex labels
        // (should be generalized, dof not always at vertices!)
        std::vector<int> loc2glb(dofs_per_elem);
        std::vector<int> dirichlet_lbs(dofs_per_elem);
        
        for (int i = 0; i < dofs_per_elem; i++) {
            loc2glb[i] = K.elem_vertices[i]->glb_idx;
            dirichlet_lbs[i] = K.elem_vertices[i]->vertex_label;    
        }

        // Loop over trial function dofs
        for (int i = 0; i < dofs_per_elem; i++) {    

            // Loop over test function dofs
            for (int j = 0; j < dofs_per_elem; j++) {

                double local_contribution = 0.;

                // Add integral over element
                for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

                    const Rn xq(qr[ipq].node);   // quadrature point in reference element

                    if (mass) {
                        psi.eval(xq, psi_vals);
                        local_contribution += mass * psi_vals[i] * psi_vals[j] * qr[ipq].weight * measure;
                    }
                                            
                    if (stiffness) {
                        psi.eval_d(K, xq, dpsi_vals);
                        for (int dm = 0; dm < dim; dm++) {
                            local_contribution += stiffness * dpsi_vals[i][dm] * dpsi_vals[j][dm] * qr[ipq].weight * measure;
                        }
                    }
                }

                if (i == j)
                    avg_diag += local_contribution / dofs_per_elem;

                // Check if either i or j is a Dirichlet dof
                int is_dirichlet = (std::find(dirichlet_lbls.begin(), dirichlet_lbls.end(), dirichlet_lbs[i]) != dirichlet_lbls.end()) || (std::find(dirichlet_lbls.begin(), dirichlet_lbls.end(), dirichlet_lbs[j]) != dirichlet_lbls.end());
                
                if (is_dirichlet && !dirichlet_lbls.empty()) {

                    local_contribution = 0.0;   // local row and column are set to zero

                    // diagonal entry is set to average of all other local diagonal entries (for better conditioning of the matrix, it could just be 1)
                    if (i == j) {

                        local_contribution = avg_diag / dofs_per_elem;
                        
                        if (!rhs.empty()) 
                            rhs[loc2glb[i]] = 0.0;
                        else {
                            std::cout << "Error! Must assemble right hand side before assembling matrix!\n";
                            exit(1);
                        }
                    }                    
                } 
                
                mat[std::make_pair(loc2glb[i], loc2glb[j])] += local_contribution;

            }

        }
    
    }

}



template <typename mesh>
template <int d, int deg>
void problem<mesh>::assemble_rhs(const QuadratureRule<d> &qr, const basis_function<d, deg> &psi, const std::function<double(const typename mesh::Rn &)> & f) {

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int ndofs = psi.ndof;             // number of dofs per element

    std::vector<double> psi_vals(ndofs);    // container for evaluations of psi

    rhs.resize(Th->nv, 0.);                 // resize rhs to number of dofs

    // Loop over all elements
    for (int k = 0; k < Th->nk; k++) {

        // Get the current element
        const elem & K = (*Th)[k];

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
void problem<mesh>::set_dirichlet(const std::function<double(const typename mesh::Rn &)> & g) {

    for (int i = 0; i < Th.nv; i++) {
        if (Th.mesh_vertices[i].vertex_label != 0) {
            rhs[i] = g(Th.mesh_vertices[i].x);

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


