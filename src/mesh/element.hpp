#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <fstream>
#include <assert.h>
#include <memory>
#include "vertex.hpp"


class Edge {

// Private variables
private:

    static const int n_vertices = 2;
    std::array<std::shared_ptr<Vertex>, n_vertices> vertices;

    int border_label;

// Constructor and methods
public:
    
    Edge() {
        for (int i = 0; i < n_vertices; i++) {
            vertices[i] = std::make_shared<Vertex>(0.);
        }
    }
    
    Edge(std::array<Vertex, n_vertices> _vertices) {
        for (int i = 0; i < n_vertices; i++) {
            vertices[i] = std::make_shared<Vertex>(_vertices[i]);
        }
    }

    Edge &set(std::vector<Vertex> v0, std::vector<int> iv, int lab) {
        for (int i = 0; i < n_vertices; ++i) 
            vertices[i] = std::make_shared<Vertex>(v0.at(iv.at(i)));

        border_label = lab;

        return *this;
    }

    const Vertex &get_vertex(unsigned int i) const {

        return *vertices.at(i);
    }

};


// Quadrilateral element structure
class Element {


// Public variables
public:
  
    static const int n_vertices = 4;    // number of vertices per element
    static const int n_edges = 4;       // number of edges per element

    int element_label;

// Private objects
private:
    std::array<std::shared_ptr<Vertex>, n_vertices> vertices;
    std::array<std::shared_ptr<Edge>, n_edges> edges;


// Constructors
public:

    Element() {
        for (int i = 0; i < n_vertices; i++) {
            vertices[i] = std::make_shared<Vertex>(0.);
        }
    }

    Element(std::array<Vertex, n_vertices> _vertices, std::array<Edge, n_edges> _edges) {
        for (int i = 0; i < n_vertices; i++) {
            vertices[i] = std::make_shared<Vertex>(_vertices[i]);
            edges[i] = std::make_shared<Edge>(_edges[i]);
        }
    }

    // Get vertex
    const Vertex &get_vertex(unsigned int i) const {
        
        return *vertices.at(i);
    }

    // Get face
    const Edge &get_edge(unsigned int i) const {

        return *edges.at(i);
    }

    Element &set(std::vector<Vertex> v0, std::vector<int> iv, int lab) {
        //assert(iv.size() == n_vertices);
        for (int i = 0; i < n_vertices; ++i) 
            vertices.at(i) = std::make_shared<Vertex>(v0.at(iv.at(i)));

        element_label = lab;

        return *this;
    }


    // Get area
    double get_area() const
    {
        return vertices[1]->x * vertices[2]->y; // base*height
    }


};

#endif // ELEMENT_HPP
