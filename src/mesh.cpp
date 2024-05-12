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


// Todo: implement Mesh2D constructor   