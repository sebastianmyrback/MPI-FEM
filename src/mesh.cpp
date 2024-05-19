#include "../include/mesh/mesh.hpp"

Mesh1D::Mesh1D(
    const double a, 
    const double b, 
    const int n) 
    {

    nverts = n + 1;
    ncells = n;
    nbe = 2;
    nsubdomains = 1;

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

    assert(false);


}


void Mesh1D::partition(const size_t n_mpi_processes) {

    // Compute local indices for data distribution
    const size_t nblocks = ncells / n_mpi_processes;      // note: integer division
    const size_t nremain = ncells % n_mpi_processes;      // remainder
    
    size_t process = 0;
    size_t cells_for_this_process = 0;

    for (int i = 0; i < ncells; ++i) 
    {
        cells[i].set_subdomain(process);

        cells_for_this_process++;
        if (cells_for_this_process >= nblocks + (process < nremain ? 1 : 0)) {
            ++process;
            cells_for_this_process = 0;
        }
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


Mesh1D::Mesh1D(
    const double a, 
    const double b, 
    const int n,
    const size_t n_mpi_processes) 
{
    nverts = n + 1;
    ncells = n;
    nbe = 2;
    nsubdomains = 1;

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

    // Partitioning
    const size_t nblocks = ncells / n_mpi_processes;      // note: integer division
    const size_t nremain = ncells % n_mpi_processes;      // remainder
    
    size_t process = 0;
    size_t cells_for_this_process = 0;

    for (int i = 0; i < ncells; ++i) 
    {
        cells[i].set_subdomain(process);

        cells_for_this_process++;
        if (cells_for_this_process >= nblocks + (process < nremain ? 1 : 0)) {
            ++process;
            cells_for_this_process = 0;
        }
    }

    nsubdomains = n_mpi_processes;

    // Distribution of DOFs
    dof_distribution.resize(nsubdomains);
    shared_dofs.clear();

    std::vector<std::unordered_set<size_t>> dof_sets(nsubdomains);

    for (int i = 0; i < ncells; ++i) 
    {
        int subdomain = cells[i].get_subdomain();

        for (int n = 0; n < cells[i].n_verts_per_cell; n++) {
            const size_t glb_idx = cell_to_vertex[cells[i].n_verts_per_cell*cells[i].get_index() + n];

            // add glb_idx to correct subdomain set if it's not already there
            dof_sets[subdomain].insert(glb_idx);
            
            shared_dofs[glb_idx].insert(subdomain);
        }
    }

    // Convert unordered_sets to vectors and sort them
    for (size_t i = 0; i < nsubdomains; ++i) {
        dof_distribution[i].assign(dof_sets[i].begin(), dof_sets[i].end());
        std::sort(dof_distribution[i].begin(), dof_distribution[i].end());
    }
}