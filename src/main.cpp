
// create a file that will be used to test the library
#include <iostream>
#include "visualization/gnuplot.hpp"
#include "mesh/mesh.hpp"
#include "mesh/element.hpp"
#include <assert.h>

int main() {
    
    // Create a mesh object
    int nx = 10, ny = 10;
    double x0 = 0.0, y0 = 0.0;
    double lx = 1.0, ly = 1.0;
    Mesh mesh(nx, ny, x0, y0, lx, ly);

    
    return 0;
}