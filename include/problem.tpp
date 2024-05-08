#include "../src/quadrature.cpp"

template <typename mesh>
template <int d, int deg>
void problem<mesh>::assemble_FEM_matrix(
    const QuadratureRule<d> &qr, 
    const basis_function<d, deg> &psi, 
    const double mass, 
    const double stiffness, 
    const dirichlet_bc<d> &bc) 
    {

    // clear the matrix
    if (!mat.empty())
        mat.clear();

    // alpha = 1 -> mass matrix
    // beta  = 1 -> stiffness matrix

    static const int dim = d;
    const int n_quad_pts = qr.n;

    const int dofs_per_elem = psi.ndof;

    std::vector<double> psi_vals(dofs_per_elem);
    std::vector<std::vector<double>> dpsi_vals(dofs_per_elem, std::vector<double>(dim, 0.));    // dofs_per_element x space dim
    
    std::map<int, double> border_dofs_values;


    // Loop over all cells in Th using the element iterator
    //for (const auto & cell: *Th)
    for (auto cell = Th->cell_begin(); cell != Th->cell_end(); ++cell)
    {
    
        double avg_diag = 0.;    // average of diagonal entries

        // Get the measure of the cell
        const double measure = cell->get_measure();

        // Create map from local to global dofs and corresponding vertex labels
        // (should be generalized, dof not always at vertices!)
        std::vector<int> loc2glb(dofs_per_elem);
        std::vector<int> labels(dofs_per_elem);
        
        for (int i = 0; i < dofs_per_elem; i++) {
            loc2glb[i] = cell->vertex(i).global_index();
            //labels[i]  = K.elem_vertices[i]->vertex_label;    

            if (bc.set_dirichlet) {

                // Add boundary value and corresponding global dof index to map if dof is on the boundary
                if (std::find(bc.lbs.begin(), bc.lbs.end(), cell->vertex(i).boundary_label()) != bc.lbs.end()) {
                    border_dofs_values.insert({cell->vertex(i).global_index(), bc.g(cell->vertex(i))});
                }

            }
        }

        // Loop over trial function dofs
        for (int i = 0; i < dofs_per_elem; i++) {    
            
            bool is_row_dof_dirichlet = false;

            
            if (bc.set_dirichlet)
                is_row_dof_dirichlet = (std::find(bc.lbs.begin(), bc.lbs.end(), cell->vertex(i).boundary_label()) != bc.lbs.end());

            // Loop over test function dofs
            for (int j = 0; j < dofs_per_elem; j++) {
                
                bool is_col_dof_dirichlet = false;

                double ak_ij = 0.;  // local matrix entry

                if (bc.set_dirichlet)
                    is_col_dof_dirichlet = (std::find(bc.lbs.begin(), bc.lbs.end(), cell->vertex(j).boundary_label()) != bc.lbs.end());

                // assemble local matrix Ak without boundary conditions 
                // subtract boundary conditions from rhs vector
                for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

                    const Rn xq(qr[ipq].node);   // quadrature point in reference element

                    if (mass) {
                        psi.eval(xq, psi_vals);
                        ak_ij += mass * psi_vals[i] * psi_vals[j] * qr[ipq].weight * measure;
                    }
                                            
                    if (stiffness) {
                        psi.eval_d(*cell, xq, dpsi_vals);

                        for (int dm = 0; dm < dim; dm++) {
                            ak_ij += stiffness * qr[ipq].weight * measure * dpsi_vals[i][dm] * dpsi_vals[j][dm];

                            if (is_col_dof_dirichlet) {
                                // subtract dirichlet column of A times boundary value from rhs
                                assert(bc.set_dirichlet);
                                rhs[loc2glb[i]] -= stiffness * qr[ipq].weight * measure * dpsi_vals[i][dm] * dpsi_vals[j][dm] * border_dofs_values.at(loc2glb[j]);
                            }
                        }
                    }
                }

                // Sum averages of diagonal entries
                if (i == j)
                    avg_diag += ak_ij / dofs_per_elem;


                // Check if either i or j is a Dirichlet dof
                bool is_dirichlet = is_row_dof_dirichlet || is_col_dof_dirichlet;
                
                if (is_dirichlet) {

                    ak_ij = 0.0;   // local row and column are set to zero

                    // diagonal entry is set to average of all other local diagonal entries (for better conditioning of the matrix, it could just be 1)
                    if (i == j) {

                        ak_ij = avg_diag;
                        
                        if (!rhs.empty()) {
                            assert(bc.set_dirichlet);
                            rhs[loc2glb[i]] = avg_diag*border_dofs_values.at(loc2glb[j]);
                        }
                        else {
                            std::cerr << "Error! Must assemble right hand side before assembling matrix!\n";
                            exit(1);
                        }
                    }                    
                } 
                
                mat[std::make_pair(loc2glb[i], loc2glb[j])] += ak_ij;

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

    rhs.resize(Th->get_nverts());                 // resize rhs to number of dofs

    // Loop over all elements
    for (auto cell = Th->cell_begin(); cell != Th->cell_end(); ++cell) {

        // Get the measure of the element
        const double measure = cell->get_measure();

        // Create map from local to global dofs
        std::vector<int> loc2glb(ndofs);
        for (int i=0; i<ndofs; i++)
            loc2glb[i] = cell->vertex(i).global_index();
        
        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts; ++ipq) {

            const Rn xq(qr[ipq].node);   // quadrature point in reference element

            Rn x;
            cell->map_to_physical(xq, x);           // map xq to quadrature point in physical element x

            psi.eval(xq, psi_vals);                     // evaluate psi at xq -> store in psi_vals
            
            // Loop over dofs
            for (int i = 0; i < ndofs; i++) {

                rhs[loc2glb[i]] += qr[ipq].weight * measure * f(x) * psi_vals[i];
            }

        }

    }

}



// template <typename mesh>
// void problem<mesh>::set_dirichlet(const std::function<double(const typename mesh::Rn &)> & g) {

//     for (int i = 0; i < Th.nv; i++) {
//         if (Th.mesh_vertices[i].vertex_label != 0) {
//             rhs[i] = g(Th.mesh_vertices[i].x);

//             // Zero out the row and column corresponding to the Dirichlet boundary
//             for (auto & [key, value] : mat) {
//                 if (key.first == i || key.second == i) {
//                     value = 0.0;
//                 }
//             }
//             mat[std::make_pair(i, i)] = 1.0;
//         }
//     }

// }


