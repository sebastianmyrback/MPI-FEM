#include <iostream>
#include "mesh.hpp"
#include <assert.h>
#include "fe.hpp"

int main() {
    
    // Create a mesh object
    int n = 10;
    double a = 0.0, b = 1.0;

    mesh1d Th(a, b, n);

    for (int i=0; i<Th.nv; i++) {
        std::cout << Th(i) << " ";
    }
    std::cout << "\n";

    std::cout << Th[0].n_vertices << std::endl;

    for (int i=0; i<Th.nk; i++) {
        //std::cout << *(Th.elements[i].elem_vertices[0]) << " ";
        std::cout << Th[i](0) << " ";
    }
    std::cout << "\n";




    //gnuplot::save(Th);

    return 0;
}