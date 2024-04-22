
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

    return 0;
}