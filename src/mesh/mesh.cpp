#include "mesh.hpp"

mesh1d::mesh1d(double a, double b, int n) {

    // n = number of vertices
    // n-1 = number of elements

    nv = n;
    nk = n-1;

    h = (b - a) / (nk);

    // Create vertices
    for (int i = 0; i < nv; i++) {
        double x = a + i * h;  // x coordinate
        mesh_vertices.push_back(vertex<1>(x));
    }

    // Create elements
    for (int i = 0; i < nk; i++) {
        std::vector<std::shared_ptr<vertex<1>>> vertices_for_element;
        vertices_for_element.push_back(std::make_shared<vertex<1>>(mesh_vertices[i]));
        vertices_for_element.push_back(std::make_shared<vertex<1>>(mesh_vertices[i + 1]));
        elements.push_back(element<1>{vertices_for_element});
    }

    return;
};