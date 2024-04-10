#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <string>

#include "element.hpp"

class Mesh {

    typedef Element::Vertex Vertex;
    typedef Element::Edge Edge;

    public:

        // Constructor
        Mesh();
        // Constructor
        Mesh(int nx, int ny, double x0, double y0, double lx, double ly);
        // Destructor
        ~Mesh();
        // Get the number of elements
        int get_n_elements();
        // Get the number of nodes
        int get_n_nodes();

        // Get the nodes
        std::vector<Vertex> get_nodes();
        // Get the elements
        std::vector<Element> get_elements();
        // Get the nodes of an element
        std::vector<Vertex> get_element_nodes(int element_id);

    private:
        
        // Number of vertices
        int n_vertices;
        // Number of elements
        int n_elements;
        // Number of border elements
        int n_be;
        // Number of nodes in the x direction
        int nx;
        // Number of nodes in the y direction
        int ny;
        // x0 value
        double x0;
        // y0 value
        double y0;
        // lx value
        double lx;
        // ly value
        double ly;
        // Nodes
        std::vector<Vertex> vertices;
        // Elements
        std::vector<Element> elements;
        // Border elements
        std::vector<Edge> borderelements;

};


#endif // MESH_HPP