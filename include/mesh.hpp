#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <array>
#include <list>
#include <memory>
#include <map>
#include <iostream>
#include <assert.h>


template<int dim>
struct Mesh;  // Forward declaration of Mesh

template<int dim>
class Point {

    double components[dim];

public:
    // Constructors
    Point() {
        for (int i = 0; i < dim; i++) {
            components[i] = 0.0;
        }
    }

    Point(const double &s) {
        for (int i = 0; i < dim; i++) {
            components[i] = s;
        }
    }

    Point(const double *s) {
        for (int i = 0; i < dim; i++) {
            components[i] = s[i];
        }
    }

    Point(const Point<dim> &p) {
        for (int i = 0; i < dim; i++) {
            components[i] = p.components[i];
        }
    }

    // Public methods

    int size() const {return dim;}


    // Operator overloading

    double& operator[](int i) {
        assert(i >= 0 && i < dim);
        return components[i];
    }

    const double& operator[](int i) const {
        assert(i >= 0 && i < dim);
        return components[i];
    }

    Point<dim> operator+(const Point<dim> &other) const {
        Point<dim> result;
        for (int i = 0; i < dim; i++) {
            result.components[i] = this->components[i] + other.components[i];
        }
        return result;
    }

    Point<dim> operator-(const Point<dim> &other) const {
        Point<dim> result;
        for (int i = 0; i < dim; i++) {
            result.components[i] = this->components[i] - other.components[i];
        }
        return result;
    }

    Point<dim> operator*(double scalar) const {
        Point<dim> result;
        for (int i = 0; i < dim; i++) {
            result.components[i] = this->components[i] * scalar;
        }
        return result;
    }

    Point<dim> operator/(double scalar) const {
        Point<dim> result;
        for (int i = 0; i < dim; i++) {
            result.components[i] = this->components[i] / scalar;
        }
        return result;
    }

    friend Point<dim> operator*(double scalar, const Point<dim> &v) {
        return v * scalar;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point<dim> &p) {
        for (int i = 0; i < dim; i++) {
            os << p.components[i] << ' ';
        }
        return os;
    }
};


template<int dim>
class Vertex : public Point<dim> {

        
    int vertex_label, glb_idx;      // label and index

public:
    int boundary_label() const {return vertex_label;}
    int global_index() const {return glb_idx;}

    Vertex() : Point<dim>() {}

    Vertex(const Point<dim> &p) : Point<dim>(p), vertex_label(0), glb_idx(0) {}

    Vertex(const Point<dim> &p, const int idx, const int label) : Point<dim>(p), glb_idx(idx), vertex_label(label) {}
};


template<int dim>
struct Cell {

    static const int n_verts_per_cell = 1 << dim;   // number of vertices in a quadrilateral cell = 2^d

    // Access the vertex with local index vertex 
    Vertex<dim> vertex(int vertex) const {
        assert(vertex >= 0 && vertex < n_verts_per_cell);
        assert(index >= 0 && index < msh->get_ncells());

        int vertex_index = msh->get_cell_to_vertex()[index*n_verts_per_cell + vertex];

        return msh->get_vertices()[vertex_index];
    }

    // Map a point from the reference element xref to the physical element X
    void map_to_physical(const Point<dim> &xref, Point<dim> &x) const {
        
        x = vertex(0) + xref[0] * (vertex(1) - vertex(0));
    };

    Cell() : msh(nullptr), index(0), measure(0.0) {assert(false);}
    Cell(Mesh<dim> *_msh, int i, double m) : msh(_msh), index(i), measure(m) {}

    int get_index() const {return index;}
    double get_measure() const {return measure;}

protected:

    Mesh<dim> *msh;

    int index;      // index of the cell in the mesh
    double measure; // measure of the cell (length, area, volume, etc.)



};


template<int dim>
struct Mesh {

protected:

    int nverts, ncells, nbe;                      // number of vertices, quadrilaterals and border elements in the mesh
    double h;                                     // typical mesh size

    std::list<Cell<dim>> cells;                     // list of cells in the mesh
    std::vector<int> cell_to_vertex;                // map index from cell to vertex
    std::vector<Vertex<dim>> vertices;                // list of vertices in the mesh

public:

    typedef Cell<dim> Element;
    typedef Point<dim> Rn;

    int get_nverts() const {return nverts;}
    int get_ncells() const {return ncells;}
    int get_nbe() const {return nbe;}
    double get_h() const {return h;}

    const std::vector<Vertex<dim>> & get_vertices() const {return vertices;}
    const std::vector<int>& get_cell_to_vertex() const {return cell_to_vertex;}


    // Iterators for iterating over the cells in the mesh
    using cell_iterator = typename std::list<Cell<dim>>::const_iterator;  

    cell_iterator cell_begin() const {
        return cells.begin();
    }

    cell_iterator cell_end() const {
        return cells.end();
    }

    using vertex_iterator = typename std::vector<Vertex<dim>>::const_iterator;

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


#endif