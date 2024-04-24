
// create a file that will be used to test the library
#include <iostream>
//#include "visualization/gnuplot.hpp"
#include "mesh.hpp"
#include <assert.h>

int main() {
    
    // Create a mesh object
    int n = 10;
    double a = 0.0, b = 1.0;

    mesh1d Th(a, b, n);

    for (int i=0; i<Th.nv; i++) {
        std::cout << Th(i) << " ";
    }
    std::cout << "\n";


    for (int i=0; i<Th.nk; i++) {
        //std::cout << *(Th.elements[i].elem_vertices[0]) << " ";
        std::cout << Th[i](0) << " ";
    }
    std::cout << "\n";




    //gnuplot::save(Th);

    return 0;
}