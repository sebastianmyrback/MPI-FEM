#include "../include/mesh/mesh.hpp"

Mesh1D::Mesh1D(
    const double a, 
    const double b, 
    const int n) 
    {

    nverts = n + 1;
    ncells = n;
    nbe = 2;

    constexpr int n_verts_per_cell = 2;

    h = (b - a) / ncells;

    vertices.reserve(nverts);
    cells.reserve(ncells);
    cell_to_vertex.assign(n_verts_per_cell * ncells, 0);
    
    vertices.push_back(Vertex<1>(Point<1>(a), 0, 1));

    // Create inner vertices
    for (int i = 1; i < nverts - 1; i++) {
        const Point<1> x(a + i*h);  
        vertices.push_back(Vertex<1>(x, i, 0));        
    }

    vertices.push_back(Vertex<1>(Point<1>(b), nverts-1, 2));

    // Create elements
    for (int i = 0; i < ncells; i++) {

        for (int j = 0; j < n_verts_per_cell; j++) {
            cell_to_vertex[i * n_verts_per_cell + j] = i + j;
        }

        cells.push_back(Cell<1>(this, i, h));    
    }


    return;

}


void Mesh1D::mesh_info(const bool detailed) 
{

    std::cout << "Number of vertices: \t \t" << nverts << std::endl;
    std::cout << "Number of cells: \t \t" << ncells << std::endl;
    std::cout << "Number of boundary elements: \t" << nbe << std::endl;
    std::cout << "Number of MPI processes: \t" << nsubdomains << std::endl;
    std::cout << "Mesh dimension: \t \t" << dim << std::endl;
    std::cout << "Mesh size: \t \t \t" << h << std::endl;
    std::cout << std::endl;
    

    if (detailed) 
    {
        for (auto cell = this->cell_begin(); cell != this->cell_end(); ++cell) 
        {
            std::cout << "Cell index " << cell->get_index() << std::endl;
            std::cout << "Cell subdomain " << cell->get_subdomain() << std::endl;
            std::cout << "Cell vertex 1 " << cell->vertex(0) << std::endl;
            std::cout << "Cell vertex 2 " << cell->vertex(1) << std::endl;
            std::cout << "Cell measure " << cell->get_measure() << std::endl;
            std::cout << std::endl;
        }
    }



}


void Mesh1D::refine_mesh() {

    // Given an existing Mesh1D object, refine the mesh by adding a new vertex in the middle of each cell
    // and creating two new cells from each old cell

    const int n_old_cells = ncells;
    const int n_old_verts = nverts;

    nverts = 2 * n_old_cells + 1;
    ncells = 2 * n_old_cells;
    nbe = 2;

    constexpr int n_verts_per_cell = 2;

    std::vector<Vertex<1>> new_vertices;
    std::vector<Cell<1>> new_cells;
    std::vector<std::size_t> new_cell_to_vertex;

    new_vertices.reserve(nverts);
    new_cells.reserve(ncells);
    new_cell_to_vertex.assign(n_verts_per_cell * ncells, 0);

    // Create new vertices
    new_vertices.push_back(vertices[0]);

    for (int i = 0; i < n_old_cells; i++) {

        const Point<1> x0 = vertices[i];
        const Point<1> x1 = vertices[i+1];
        const Point<1> xmid = 0.5 * (x0 + x1);

        new_vertices.push_back(Vertex<1>(xmid, i, 0));
        new_vertices.push_back(vertices[i+1]);

    }

    // Create new cells
    for (int i = 0; i < n_old_cells; i++) {

        const std::size_t v0 = i;
        const std::size_t v1 = i + 1;
        const std::size_t v2 = i + 2;

        const std::size_t v3 = i + 2;
        const std::size_t v4 = i + 3;
        const std::size_t v5 = i + 4;

        new_cell_to_vertex[2*i] = v0;
        new_cell_to_vertex[2*i + 1] = v2;

        new_cells.push_back(Cell<1>(this, 2*i, h/2));
        new_cell_to_vertex[2*i + 2] = v3;
        new_cell_to_vertex[2*i + 3] = v5;

        new_cells.push_back(Cell<1>(this, 2*i + 1, h/2));

    }

    vertices = new_vertices;
    cells = new_cells;
    cell_to_vertex = new_cell_to_vertex;



    return;

}


void Mesh1D::partition(const size_t n_mpi_processes) {

    // Compute local indices for data distribution
    const size_t nblocks = ncells / n_mpi_processes;      // note: integer division
    const size_t nremain = ncells % n_mpi_processes;      // remainder
    
    size_t process = 0;
    size_t cells_for_this_process = 0;

    dof_distribution.resize(nsubdomains);
    shared_dofs.clear();

    //for (auto cell = this->cell_begin(); cell != this->cell_end(); ++cell)
    for (int i = 0; i < ncells; ++i) 
    {
        cells[i].set_subdomain(process);

        cells_for_this_process++;
        if (cells_for_this_process >= nblocks + (process < nremain ? 1 : 0)) {
            ++process;
            cells_for_this_process = 0;
        }

        // for (int n = 0; n < cells[i].n_verts_per_cell; n++) {
        //     const size_t glb_idx = cell_to_vertex[cells[i].n_verts_per_cell*cells[i].get_index() + n];

        //     // add glb_idx to correct subdomain list if it's not already there
        //     if (std::find(dof_distribution[cells[i].get_subdomain()].begin(), dof_distribution[cells[i].get_subdomain()].end(), glb_idx) == dof_distribution[cells[i].get_subdomain()].end())
        //         dof_distribution[cells[i].get_subdomain()].push_back(glb_idx);
            
        //     shared_dofs[glb_idx].insert(process);
        // }


    }

    nsubdomains = n_mpi_processes;

    return;
}


void Mesh1D::distribute_dofs()
{
    dof_distribution.resize(nsubdomains);
    shared_dofs.clear();

    //for (auto cell = this->cell_begin(); cell != this->cell_end(); ++cell)
    for (int i = 0; i < ncells; ++i) 
    {
        int subdomain = cells[i].get_subdomain();

        for (int n = 0; n < cells[i].n_verts_per_cell; n++) {
            const size_t glb_idx = cell_to_vertex[cells[i].n_verts_per_cell*cells[i].get_index() + n];

            // add glb_idx to correct subdomain list if it's not already there
            if (std::find(dof_distribution[cells[i].get_subdomain()].begin(), dof_distribution[cells[i].get_subdomain()].end(), glb_idx) == dof_distribution[cells[i].get_subdomain()].end())
                dof_distribution[cells[i].get_subdomain()].push_back(glb_idx);
            
            shared_dofs[glb_idx].insert(subdomain);
        }

    }

}

// const std::vector<std::vector<size_t>> &Mesh1D::get_distribution() {

//     // Return a vector of vectors containing the global indices
//     // of the local dofs on each process

//     std::vector<std::vector<size_t>> dof_distribution(nsubdomains);

//     for (auto cell = this->cell_begin(); cell != this->cell_end(); ++cell) {

//         for (int n = 0; n < cell->n_verts_per_cell; n++) {
//             const size_t glb_idx = cell_to_vertex[cell->n_verts_per_cell*cell->get_index() + n];

//             // add glb_idx to correct subdomain list if it's not already there
//             if (std::find(dof_distribution[cell->get_subdomain()].begin(), dof_distribution[cell->get_subdomain()].end(), glb_idx) == dof_distribution[cell->get_subdomain()].end())
//                 dof_distribution[cell->get_subdomain()].push_back(glb_idx);
            
//         }
//     }

//     return dof_distribution;

// }


// const std::map<int, std::set<int>> &get_shared_dofs() {

//     // Return a map of shared dofs for parallelization

//     std::map<int, std::set<int>> shared_dofs;

//     for (auto cell = this->cell_begin(); cell != this->cell_end(); ++cell) {
//         int subdomain = cell->get_subdomain();
//         for (size_t i = 0; i < dofs_per_cell; ++i) {
//             int dof = cell->vertex(i).global_index();
//             shared_dofs[dof].insert(subdomain);
//         }
//     }

//     return shared_dofs;

// }