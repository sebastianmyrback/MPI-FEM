#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <assert.h>

template<int dim>
struct Mesh;  // Forward declaration of Mesh


template<int dim>
class Point {
    std::array<double, dim> components;

public:

    Point() : components{} {}

    explicit Point(const double &s) : components{s} {}

    explicit Point(const std::array<double, dim> &s) : components{s} {}

    Point(const Point<dim> &p) : components{p.components} {}

    std::size_t size() const noexcept { return dim; }

    double& operator[](std::size_t i) {
        assert(i < dim);
        return components[i];
    }

    const double& operator[](std::size_t i) const {
        assert(i < dim);
        return components[i];
    }

    Point operator+(const Point &other) const;
    Point operator-(const Point &other) const;
    Point operator*(double scalar) const;
    Point operator/(double scalar) const;

    // friend methods must be defined in the class block
    // since the class is a template and hence only
    // is created at run-time.

    friend Point operator*(double scalar, const Point &v) 
    {
        return v * scalar;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &p) 
    {
        for (std::size_t i = 0; i < dim; i++) {
            os << p.components[i] << ' ';
        }
        return os;
    }

};



template<int dim>
class Vertex : public Point<dim> {

        
    std::size_t vertex_label, glb_idx;      // label and index

public:
    std::size_t boundary_label() const {return vertex_label;}
    std::size_t global_index() const {return glb_idx;}

    Vertex() : Point<dim>() {}

    Vertex(const Point<dim> &p) : Point<dim>(p), vertex_label(0), glb_idx(0) {}

    Vertex(const Point<dim> &p, const std::size_t idx, const std::size_t label) : Point<dim>(p), glb_idx(idx), vertex_label(label) {}
};


template<int dim>
struct Cell {

    static constexpr std::size_t n_verts_per_cell = 1 << dim;   // number of vertices in a quadrilateral cell = 2^d

    // Access the vertex with local index vertex 
    Vertex<dim> vertex(std::size_t vertex) const {
        assert(vertex >= 0 && vertex < n_verts_per_cell);
        assert(index >= 0 && index < msh->get_ncells());

        std::size_t vertex_index = msh->get_cell_to_vertex()[index*n_verts_per_cell + vertex];

        return msh->get_vertices()[vertex_index];
    }

    // Map a point from the reference element xref to the physical element X
    void map_to_physical(const Point<dim> &xref, Point<dim> &x) const {
        
        x = vertex(0) + xref[0] * (vertex(1) - vertex(0));
    };

    Cell() : msh(nullptr), index(0), measure(0.0) {assert(false);}
    Cell(Mesh<dim> *_msh, std::size_t i, double m) : msh(_msh), index(i), measure(m) {}

    std::size_t get_index() const {return index;}
    double get_measure() const {return measure;}

protected:

    const Mesh<dim> *msh;

    std::size_t index;      // index of the cell in the mesh
    double measure; // measure of the cell (length, area, volume, etc.)

};


template<int d>
struct Mesh {

protected:

    std::size_t nverts, ncells, nbe;               // number of vertices, quadrilaterals and border elements in the mesh
    double h;                              // typical mesh size

    std::vector<Cell<d>> cells;            // list of cells in the mesh
    std::vector<Vertex<d>> vertices;     // list of vertices in the mesh
    std::vector<std::size_t> cell_to_vertex;       // map index from cell to vertex

public:

    const std::size_t get_nverts() const {return nverts;}
    const std::size_t get_ncells() const {return ncells;}
    const std::size_t get_nbe() const {return nbe;}
    const std::size_t get_dim() const {return d;}
    const double get_h() const {return h;}


    const std::vector<Vertex<d>> & get_vertices() const {return vertices;}
    const std::vector<std::size_t>& get_cell_to_vertex() const {return cell_to_vertex;}


    // Iterators for iterating over the cells in the mesh
    using cell_iterator = typename std::vector<Cell<d>>::const_iterator;  

    cell_iterator cell_begin() const {
        return cells.begin();
    }

    cell_iterator cell_end() const {
        return cells.end();
    }

    using vertex_iterator = typename std::vector<Vertex<d>>::const_iterator;

    vertex_iterator vertex_begin() const {
        return vertices.begin();
    }

    vertex_iterator vertex_end() const {
        return vertices.end();
    }

};



struct Mesh1D : public Mesh<1> {

    Mesh1D(const double a, const double b, const int n);
};

struct Mesh2D : public Mesh<2> {

    Mesh2D(const double x0, const double y0, const double xend, const double yend, const int nx, const int ny);

};

#include "mesh.tpp"

