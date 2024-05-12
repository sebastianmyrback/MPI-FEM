#include "mesh.hpp"

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
    for (int i = 0; i < ncells; ++i) 
    {
        cells[i].set_subdomain(process);

        if ((i + 1) % nblocks == 0 && (nremain == 0 || process >= nremain)) {
            ++process;
        }
    }

    return;

}

