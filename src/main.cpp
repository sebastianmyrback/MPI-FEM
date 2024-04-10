
// create a file that will be used to test the library
#include <iostream>
#include "mesh/mesh.hpp"
#include "mesh/element.hpp"

int main() {
    
    // Create a mesh object
    int nx, ny = 10;
    double x0, y0 = 0.0;
    double lx, ly = 1.0;
    //Mesh mesh(nx, ny, x0, y0, lx, ly);

    typedef Element::Vertex Vertex;
    typedef Element::Edge Edge;

    // Create a vertex object
    std::array<Element::Vertex, 4> vertices = {Vertex(0.0, 0.0), Vertex(1.0, 0.0), Vertex(1.0, 1.0), Vertex(0.0, 1.0)};
    std::array<Edge, 4> edges = {Edge({vertices[0], vertices[1]}), Edge({vertices[1], vertices[2]}), Edge({vertices[2], vertices[3]}), Edge({vertices[3], vertices[0]})};

    Element element(vertices, edges);

    std::cout << edges[0][0] << " " << edges[0][1] << std::endl;
    std::cout << vertices[1] << " " << vertices[2] << std::endl;
    std::cout << vertices[2] << " " << vertices[3] << std::endl;

    
    std::cout << "Hello, World!" << std::endl;
    return 0;
}