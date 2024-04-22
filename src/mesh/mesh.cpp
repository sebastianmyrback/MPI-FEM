// Implement the mesh class

#include <assert.h>
#include "mesh.hpp"
#include "element.hpp"
#include <iostream>

// Constructor
Mesh::Mesh(){}

// Constructor
Mesh::Mesh(int _nx, int _ny, double _x0, double _y0, double _lx, double _ly) : nx(_nx), ny(_ny), x0(_x0), y0(_y0), lx(_lx), ly(_ly) {

    const std::array<int, 4> indQ = {0, 1, 3, 2};

    n_vertices = nx * ny;
    n_elements = (nx - 1) * (ny - 1);
    n_be = 2 * ((nx - 1) + (ny - 1));
    const double hx = lx / (nx - 1);
    const double hy = ly / (ny - 1);

    vertices.resize(n_vertices);
    elements.resize(n_elements);
    borderelements.resize(n_be);

    std::vector<int> idc_vert_mesh(4);  // vertex indices in the order of iterating through the whole mesh
    std::vector<int> indc_vert_elem(4); // vertex indices in the order each vertex is stored in an element

    int jt = 0;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {

            int id = 0;
            for (int jj = j; jj < j + 2; ++jj) {
                for (int ii = i; ii < i + 2; ++ii) {

                    int ivl  = ii + jj * nx; // index
                    idc_vert_mesh[id++] = ivl;

                    R2 P(x0 + ii*hx, y0 + jj*hy);
                    vertices.at(ivl) = Vertex(P);
                    
                }
            }
            
            for (int e = 0; e < 4; ++e) {
                indc_vert_elem[e] = idc_vert_mesh[indQ[e]];
            }
            
            
            elements[jt++].set(vertices, indc_vert_elem, 0);

        }

    }

    // create the for borders
    int lab, k = 0;
    for (int i = 0; i < nx - 1; ++i) {
        indc_vert_elem[0] = i;
        indc_vert_elem[1] = i + 1;
        lab = 1;
        for (int j = 0; j < 2; ++j)
            (vertices.at(indc_vert_elem[j])).vertex_label = std::max(vertices[indc_vert_elem[j]].vertex_label, lab);
        borderelements[k++].set(vertices, indc_vert_elem, lab);
    }
    for (int i = 0; i < ny - 1; ++i) {
        indc_vert_elem[0] = (i + 1) * nx - 1;
        indc_vert_elem[1] = indc_vert_elem[0] + nx;
        lab     = 2;
        for (int j = 0; j < 2; ++j)
            vertices[indc_vert_elem[j]].vertex_label = std::max(vertices[indc_vert_elem[j]].vertex_label, lab);
        borderelements[k++].set(vertices, indc_vert_elem, lab);
    }
    for (int i = 0; i < nx - 1; ++i) {
        indc_vert_elem[0] = i + nx * (ny - 1);
        indc_vert_elem[1] = indc_vert_elem[0] + 1;
        lab     = 3;
        for (int j = 0; j < 2; ++j)
            vertices[indc_vert_elem[j]].vertex_label = std::max(vertices[indc_vert_elem[j]].vertex_label, lab);
        borderelements[k++].set(vertices, indc_vert_elem, lab);
    }
    for (int i = 0; i < ny - 1; ++i) {
        indc_vert_elem[0] = i * nx;
        indc_vert_elem[1] = indc_vert_elem[0] + nx;
        lab     = 4;
        for (int j = 0; j < 2; ++j)
            vertices[indc_vert_elem[j]].vertex_label = std::max(vertices[indc_vert_elem[j]].vertex_label, lab);
        borderelements[k++].set(vertices, indc_vert_elem, lab);
    }
    
}






