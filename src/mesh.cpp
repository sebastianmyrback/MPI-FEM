#include "mesh.hpp"

#include <iostream>

mesh1d::mesh1d(double a, double b, int n) {

    // n = number of vertices
    // n-1 = number of elements

    this->nv = n+1;
    nk = n;

    h = (b - a) / (nk);

    this->mesh_vertices.resize(nv);
    this->elements.resize(nk);

    // Create vertices
    for (int i = 0; i < nv; i++) {
        const Rn x(a + i*h);  // x coordinate
        (this->mesh_vertices)[i].x = x;
        (this->mesh_vertices)[i].idx = i;
        (this->mesh_vertices)[i].vertex_label = 0;
    }

    // Create elements
    for (int i = 0; i < nk; i++) {
        std::vector<std::shared_ptr<vert>> vertices_for_element;
        vertices_for_element.push_back(std::make_shared<vert>(mesh_vertices[i]));
        vertices_for_element.push_back(std::make_shared<vert>(mesh_vertices[i + 1]));
        (this->elements)[i].elem_vertices = vertices_for_element;
        (this->elements)[i].measure = h;
    }

    // Mark the boundary vertices
    (this->mesh_vertices)[0].vertex_label = 1;
    (this->mesh_vertices)[nv-1].vertex_label = 2;

    return;
};

