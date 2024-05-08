#include "mesh.hpp"

Mesh1D::Mesh1D(
    const double a, 
    const double b, 
    const int n) 
    {

    nverts = n + 1;
    ncells = n;
    nbe = 2;

    h = (b - a) / ncells;

    vertices.resize(nverts);
    cell_to_vertex.resize(n_verts_per_cell * ncells);
    border_dofs.reserve(nbe);

    // Create inner vertices
    for (int i = 1; i < nverts - 1; i++) {
        const Rd<1> x(a + i*h);  
        vertices[i] = Vertex<1>(x, i, 0);
        
    }

    // Mark the boundary vertices
    vertices[0] = Vertex<1>(a, 0, 1);
    vertices[nverts-1] = Vertex<1>(b, nverts-1, 2);

    border_dofs.push_back(0);
    border_dofs.push_back(nverts-1);

    // Create elements
    int l = 0;
    for (int i = 0; i < ncells; i++) {

        for (int j = 0; j < n_verts_per_cell; j++) {
            cell_to_vertex[i * n_verts_per_cell + j] = i + j;
        }

        cells.push_back(Cell<1>(std::make_shared<Mesh1D>(*this), i, h));      
    
    }


    return;

}
