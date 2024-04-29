
// create a file that will be used to test the library
#include <iostream>
#include "mesh/mesh.hpp"
#include <assert.h>

int main() {
    
    Mesh mesh(3, 3, 0, 0, 2, 2);

    assert(mesh.get_n_elements() == 4);
    assert(mesh.get_n_vertices() == 9);
    assert(mesh.get_n_be() == 8);

    const Element &T(mesh.get_element(0)); // first element and the unit square
    assert(T.get_vertex(0)[0] == 0.);
    assert(T.get_vertex(0)[1] == 0.);
    assert(T.get_vertex(1)[0] == 1.);
    assert(T.get_vertex(1)[1] == 0.);
    assert(T.get_vertex(2)[0] == 1.);
    assert(T.get_vertex(2)[1] == 1.);
    assert(T.get_vertex(3)[0] == 0.);
    assert(T.get_vertex(3)[1] == 1.);

    const Vertex &V(mesh.get_vertex(4)); // vertex with index 4
    assert(V[0] == 1.);
    assert(V[1] == 1.);

    const Edge &BE(mesh.get_border_element(2)); // border element with index 2
    assert(BE.get_vertex(0)[0] == 2.);
    assert(BE.get_vertex(0)[1] == 0.);
    assert(BE.get_vertex(1)[0] == 2.);
    assert(BE.get_vertex(1)[1] == 1.);


        // Create a mesh object
    int n = 10;     // number of elements
    double a = 0.0, b = 1.0;

    mesh1d Th(a, b, n);

    for (int i=0; i<Th.nv; i++) {
        std::cout << Th(i).x << " ";
    }
    std::cout << "\n";

    std::cout << Th[0].n_vertices << std::endl;

    for (int i=0; i<Th.nk; i++) {
        //std::cout << *(Th.elements[i].elem_vertices[0]) << " ";
        std::cout << (Th[i](0)).x << " ";
    }
    std::cout << "\n";


    return 0;
}