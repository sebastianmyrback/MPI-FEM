#include "mesh.hpp"

#include <iostream>

mesh1d::mesh1d(double a, double b, int n) {

    // n = number of vertices
    // n-1 = number of elements

    nv = n;
    nk = n-1;

    h = (b - a) / (nk);

    // Create vertices
    for (int i = 0; i < nv; i++) {
        double x = a + i * h;  // x coordinate
        mesh_vertices.push_back(vert(x));
    }

    // Create elements
    for (int i = 0; i < nk; i++) {
        std::vector<std::shared_ptr<vert>> vertices_for_element;
        vertices_for_element.push_back(std::make_shared<vert>(mesh_vertices[i]));
        vertices_for_element.push_back(std::make_shared<vert>(mesh_vertices[i + 1]));
        elements.push_back(elem{vertices_for_element});
    }

    return;
};

