#include "mesh.hpp"

Mesh1D::Mesh1D(
    const double a, 
    const double b, 
    const int n) 
    {

    nverts = n + 1;
    ncells = n;
    nbe = 2;

    const int n_verts_per_cell = Cell<1>::n_verts_per_cell;

    h = (b - a) / ncells;

    vertices.resize(nverts);
    cell_to_vertex.resize(n_verts_per_cell * ncells);

    // Create inner vertices
    for (int i = 1; i < nverts - 1; i++) {
        const Point<1> x(a + i*h);  
        vertices[i] = Vertex<1>(x, i, 0);
        
    }

    // Mark the boundary vertices
    vertices[0] = Vertex<1>(a, 0, 1);
    vertices[nverts-1] = Vertex<1>(b, nverts-1, 2);


    // Create elements
    int l = 0;
    for (int i = 0; i < ncells; i++) {

        for (int j = 0; j < n_verts_per_cell; j++) {
            cell_to_vertex[i * n_verts_per_cell + j] = i + j;
        }

        //cells.push_back(Cell<1>(std::make_shared<Mesh1D>(*this), i, h));      
        cells.push_back(Cell<1>(this, i, h));    
    }


    return;

}

// Implement     Mesh2D(const double x0, const double y0, const double xend, const double yend, const int nx, const int ny);

Mesh2D::Mesh2D(
    const double x0, 
    const double y0, 
    const double xend, 
    const double yend, 
    const int nx, 
    const int ny) 
    {

    nverts = (nx + 1) * (ny + 1);
    ncells = nx * ny;
    nbe = 2 * (nx + ny);

    const int n_verts_per_cell = Cell<2>::n_verts_per_cell;

    h = (xend - x0) / nx;

    vertices.resize(nverts);
    cell_to_vertex.resize(n_verts_per_cell * ncells);

    // Create inner vertices
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny; j++) {
            double x_[2] = {x0 + i*h, y0 + j*h};
            const Point<2> x(x_);  
            vertices[i * (ny + 1) + j] = Vertex<2>(x, i * (ny + 1) + j, 0);
        }
    }

    // Mark the boundary vertices
    for (int i = 0; i < nx + 1; i++) {
        vertices[i * (ny + 1)] = Vertex<2>(x0 + i*h, i * (ny + 1), 1);
        vertices[i * (ny + 1) + ny] = Vertex<2>(x0 + i*h, i * (ny + 1) + ny, 3);
    }

    for (int j = 0; j < ny + 1; j++) {
        vertices[j] = Vertex<2>(y0 + j*h, j, 4);
        vertices[nx * (ny + 1) + j] = Vertex<2>(y0 + j*h, nx * (ny + 1) + j, 2);
    }

    // Create elements
    int l = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {

            for (int k = 0; k < n_verts_per_cell; k++) {
                cell_to_vertex[l * n_verts_per_cell + k] = i * (ny + 1) + j + k;
            }

            //cells.push_back(Cell<2>(std::make_shared<Mesh2D>(*this), l, h));
            cells.push_back(Cell<2>(this, l, h));
            l++;
        }
    }

    return;

}

