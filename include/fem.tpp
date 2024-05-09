template<typename mesh>
template <int dim, int degree>
void FEM<mesh>::compute_stiffness_on_cell(
    const typename mesh::cell_iterator &cell,
    const std::vector<int> &loc2glb,
    const QuadratureRule<dim> &qr,
    const BasisFunction<dim, degree> &psi,
    std::vector<std::vector<double>> &Ak) 
{

    // This function computes the local stiffness matrix Ak on the current cell

    const int n_quad_pts = qr.n;
    const int dofs_per_cell = psi.ndof;

    Ak.assign(dofs_per_cell, std::vector<double>(dofs_per_cell, 0.0));

    // Holder for evaluations of the gradient of psi
    std::vector<std::vector<double>> dpsi_vals(dofs_per_cell, std::vector<double>(dim, 0.));    // dofs_per_cell x space dim

    const double measure = cell->get_measure();

    // Loop over trial function dofs
    for (int i = 0; i < dofs_per_cell; ++i) 
    {        
        // Loop over test function dofs
        for (int j = 0; j < dofs_per_cell; ++j) 
        {    
            // Integrate shape functions of current dofs over the cell 
            for (int ipq = 0; ipq < n_quad_pts; ++ipq) 
            {
                const Point<dim> xq(qr[ipq].node);   // quadrature point in reference element
                psi.eval_d(*cell, xq, dpsi_vals);    // evaluate gradient of psi at xq -> store in dpsi_vals

                for (int dm = 0; dm < dim; dm++)     // loop over space dimensions
                { 
                    Ak[i][j] += qr[ipq].weight * measure * dpsi_vals[i][dm] * dpsi_vals[j][dm];
                }       
            }
        }
    }
}


template<typename mesh>
template <int dim, int degree>
void FEM<mesh>::compute_rhs_on_cell(
    const typename mesh::cell_iterator &cell,
    const std::vector<int> &loc2glb,
    const QuadratureRule<dim> &qr,
    const BasisFunction<dim, degree> &psi,
    const double f(const Point<dim> &),
    std::vector<double> &fk)
{

    // This function computes the local rhs vector fk on the current cell

    const int n_quad_pts = qr.n;
    const int dofs_per_cell = psi.ndof;

    fk.assign(dofs_per_cell, 0.0);

    std::vector<double> psi_vals(dofs_per_cell);    // container for evaluations of psi

    // Get the measure of the element
    const double measure = cell->get_measure();

    // Loop over dofs
    for (int i = 0; i < dofs_per_cell; i++) 
    {
    
        // Loop over quadrature points
        for (int ipq = 0; ipq < n_quad_pts; ++ipq) 
        {

            Point<dim> x;
            const Point<dim> xq(qr[ipq].node); // quadrature point in reference element
            cell->map_to_physical(xq, x);      // map xq to quadrature point in physical element x
            psi.eval(xq, psi_vals);            // evaluate psi at xq -> store in psi_vals

            fk[i] += qr[ipq].weight * measure * f(x) * psi_vals[i];
        }

    }

}

template<typename mesh>
template <int dim>
void FEM<mesh>::get_boundary_data(
    const typename mesh::cell_iterator &cell,
    const std::vector<int> &loc2glb,
    const DirichletBC<dim> &bc,
    std::unordered_map<int, double> &boundary_data)
{
    
    // This function checks if the current cell
    // contains any boundary dofs and adds the corresponding
    // global dof indices and boundary values to the boundary_data map.

    for (int i = 0; i < loc2glb.size(); i++) 
    {
        
        if (bc.set_dirichlet) 
        {

            // Add data if dof is a boundary dof
            if (std::find(bc.lbs.begin(), bc.lbs.end(), cell->vertex(i).boundary_label()) != bc.lbs.end())
                boundary_data.emplace(cell->vertex(i).global_index(), bc.g(cell->vertex(i)));
            
        }
    }

}

template<typename mesh>
void FEM<mesh>::distribute_local_to_global(
    const typename mesh::cell_iterator &cell,
    const std::vector<int> &loc2glb,
    std::vector<std::vector<double>> &Ak,
    std::vector<double> &fk,
    const std::unordered_map<int, double> &boundary_data)
{

    // This function assembles the local stiffness matrix Ak and local rhs vector
    // fk into the global matrix mat and global rhs vector rhs.
    // Moreover, it applies strong Dirichlet boundary conditions to the system
    // if the user has specified so in the DirichletBC object.

    if (loc2glb.size() != Ak.size() || loc2glb.size() != fk.size())
        throw std::invalid_argument("Mismatch in sizes of loc2glb, Ak, and fk");


    for (int i = 0; i < loc2glb.size(); i++)
    {

        // Add local rhs vector to global rhs vector
        rhs[loc2glb[i]] += fk[i];

        
        // Check if dof i is a Dirichlet dof
        double avg_diag = 0.; 
        bool is_row_dof_dirichlet = false;
        is_row_dof_dirichlet = (boundary_data.find(loc2glb[i]) != boundary_data.end()); 
    

        for (int j = 0; j < loc2glb.size(); j++)
        {
        
            avg_diag += Ak[j][j] / loc2glb.size();

            // Check if dof j is a Dirichlet dof
            bool is_col_dof_dirichlet = false;
            is_col_dof_dirichlet = (boundary_data.find(loc2glb[j]) != boundary_data.end());

            // If i or j is a Dirichlet dof:
            // set Dirichlet condition in local matrix and
            // change the corresponding entry in the rhs vector
            if (is_row_dof_dirichlet || is_col_dof_dirichlet) 
            {
                // Subtract dirichlet column of A times boundary value from rhs
                if (is_col_dof_dirichlet) 
                {
                    rhs[loc2glb[i]] -= Ak[i][j] * boundary_data.at(loc2glb[j]);
                }

                // Set local row and column to zero
                Ak[i][j] = 0.0;

                // Set diagonal entry to average of all other local diagonal entries
                // modify rhs accordingly
                if (i == j) 
                {
                    Ak[i][j] = avg_diag;
                    rhs[loc2glb[i]] = avg_diag * boundary_data.at(loc2glb[j]);
                }
            }

            // Add local matrix entry to global matrix
            mat[std::make_pair(loc2glb[i], loc2glb[j])] += Ak[i][j];
            
        }
    }
}

template <typename mesh>
template <int dim, int deg>
void FEM<mesh>::assemble_stiffness_system(
    const QuadratureRule<dim> &qr, 
    const BasisFunction<dim, deg> &psi, 
    const double f(const Point<dim> &),
    const DirichletBC<dim> &bc) 
{
    
    if (dim != Th->get_dim())
        throw std::invalid_argument("Mismatch in dimensions");

    if (!mat.empty())
        throw std::runtime_error("Matrix is not empty");

    
    rhs.assign(Th->get_nverts(), 0);

    const int dofs_per_cell = psi.ndof;

    // Loop over all cells in the mesh
    for (auto cell = Th->cell_begin(); cell != Th->cell_end(); ++cell) 
    {

        std::vector<int> loc2glb(dofs_per_cell);    
        for (int i = 0; i < dofs_per_cell; i++) 
            loc2glb[i] = cell->vertex(i).global_index();

        std::vector<std::vector<double>> Ak;    
        std::vector<double> fk;
        

        std::unordered_map<int, double> boundary_data;

        compute_stiffness_on_cell(cell, loc2glb, qr, psi, Ak);      
        compute_rhs_on_cell(cell, loc2glb, qr, psi, f, fk);
        get_boundary_data(cell, loc2glb, bc, boundary_data);
        distribute_local_to_global(cell, loc2glb, Ak, fk, boundary_data);

    }

}