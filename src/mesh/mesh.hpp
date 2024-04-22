#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <string>

#include "element.hpp"

class Mesh {


    public:

        Mesh();
        Mesh(int nx, int ny, double x0, double y0, double lx, double ly);
        ~Mesh() {};

        int get_n_elements() const {return n_elements;}
        int get_n_vertices() const {return n_vertices;}
        int get_n_be() const {return n_be; }

        std::vector<Vertex> get_nodes();
        std::vector<Element> get_elements();
        std::vector<Vertex> get_element_nodes(int element_id);

        const Vertex &get_vertex(int i) const {return vertices.at(i);}
        const Element &get_element(int i) const {return elements.at(i);}
        const Edge &get_border_element(int i) const {return borderelements.at(i);}
        
    private:
        
        int n_vertices = 0;
        int n_elements = 0;
        int n_be = 0;
        
        int nx, ny;
        double x0, y0;
        double lx, ly;
        
        std::vector<Vertex> vertices;
        std::vector<Element> elements;
        std::vector<Edge> borderelements;


};


#endif // MESH_HPP