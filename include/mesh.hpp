#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <array>
#include <list>
#include <memory>
#include <iostream>
#include <assert.h>


template<int d>
struct Mesh;  // Forward declaration of Mesh


template<int d>
struct Rd {
    std::vector<double> x;

    // Default constructor
    Rd() : x(d, 0.0) {}

    Rd(const double &s) : x(d, s) {}

    // Constructor that takes a std::vector as an argument
    Rd(const std::vector<double>& v) {
        if (v.size() != d) {
            throw std::invalid_argument("Vector must have exactly " + std::to_string(d) + " elements");
        }
        x = v;
    }

    // Constructor that takes a std::array as an argument
    Rd(const std::array<double, d>& a) {
        x.assign(a.begin(), a.end());
    }

    int size() const {return x.size();}

    const double & operator[](int i) const {return x.at(i); }


    // Operator overloading for arithmetic operators
    Rd operator+(const Rd &other) const {
        Rd result;
        for (int i = 0; i < d; i++) {
            result.x[i] = this->x[i] + other.x[i];
        }
        return result;
    }

    Rd operator-(const Rd &other) const {
        Rd result;
        for (int i = 0; i < d; i++) {
            result.x[i] = this->x[i] - other.x[i];
        }
        return result;
    }

    Rd& operator+=(const Rd &other) {
        for (int i = 0; i < d; i++) {
            this->x[i] += other.x[i];
        }
        return *this;
    }

    Rd& operator-=(const Rd &other) {
        for (int i = 0; i < d; i++) {
            this->x[i] -= other.x[i];
        }
        return *this;
    }

    Rd operator*(double scalar) const {
        Rd result;
        for (int i = 0; i < d; i++) {
            result.x[i] = this->x[i] * scalar;
        }
        return result;
    }

    Rd operator/(double scalar) const {
        Rd result;
        for (int i = 0; i < d; i++) {
            result.x[i] = this->x[i] / scalar;
        }
        return result;
    }

    // Operator overloading for multiplication with a scalar
    friend Rd operator*(double scalar, const Rd &v) {
        return v * scalar;
    }

    // Dot product
    double dot(const Rd &other) const {
        double result = 0.0;
        for (int i = 0; i < d; i++) {
            result += this->x[i] * other.x[i];
        }
        return result;
    }

    friend std::ostream &operator<<(std::ostream &os, const Rd<d> &rd) {
        for (const auto &val : rd.x) {
            os << val << ' ';
        }
        return os;
    }

};


template<int d>
struct Vertex : public Rd<d> {

        
    int vertex_label, glb_idx;      // label and index

    int boundary_label() const {return vertex_label;}
    int global_index() const {return glb_idx;}

    Vertex() : Rd<d>() {}

    Vertex(const Rd<d> &r) : Rd<d>(r), vertex_label(0), glb_idx(0) {}

    Vertex(const Rd<d> &r, const int idx, const int label) : Rd<d>(r), glb_idx(idx), vertex_label(label) {}
};

// Struct representing a quadrilateral element in dimension d 
template<int d>
struct Quad {

    // Access the vertex with local index vertex 
    Vertex<d> vertex(int vertex) const {
        assert(0 < vertex < n_verts_per_quad);
        assert(0 < index < msh->get_nquads());
        int vertex_index = msh->get_quad_to_vertex()[index*n_verts_per_quad + vertex];

        return msh->get_vertices()[vertex_index];
    }

    // Map a point from the reference element xref to the physical element X
    void map_to_physical(const Rd<d> &xref, Rd<d> &x) const {
        
        x = vertex(0) + xref[0] * (vertex(1) - vertex(0));
    };

    Quad() : msh(nullptr), index(0), measure(0.0) {}
    Quad(std::shared_ptr<Mesh<d>> _msh) : msh(_msh), index(0), measure(0.0) {}
    Quad(std::shared_ptr<Mesh<d>> _msh, int i) : msh(_msh), index(i), measure(0.0) {}
    Quad(std::shared_ptr<Mesh<d>> _msh, int i, double m) : msh(_msh), index(i), measure(m) {}

    int get_index() const {return index;}
    double get_measure() const {return measure;}

protected:
    
    std::shared_ptr<Mesh<d>> msh;
    int index;      // index of the quad in the mesh
    double measure; // measure of the quad (length, area, volume, etc.)

    static const int n_verts_per_quad = 1 << d;   // number of vertices in a quadrilateral quad = 2^d

};

// Struct representing a quadrilateral mesh in dimension d
template<int d>
struct Mesh {

protected:

    std::list<Quad<d>> quads;                     // list of quadrialterals
    std::vector<int> border_dofs;                 // list of global border dofs
    std::vector<int> quad_to_vertex;              // map index from quad to vertex
    std::vector<Vertex<d>> vertices;              // list of vertices in the mesh

    static const int n_verts_per_quad = 1 << d;   // number of vertices in a quadrilateral quad = 2^d

    int nverts, nquads, nbe;                      // number of vertices, quadrilaterals and border elements in the mesh
    double h;                                     // typical mesh size

public:

    typedef Quad<d> Element;
    typedef Rd<d> Rn;

    int get_nverts() const {return nverts;}
    int get_nquads() const {return nquads;}
    int get_nbe() const {return nbe;}
    double get_h() const {return h;}

    std::vector<int> get_border_dofs() const {return border_dofs;}
    std::vector<Vertex<d>> get_vertices() const {return vertices;}
    std::vector<int> get_quad_to_vertex() const {return quad_to_vertex;}


    // /**
    //  * @brief Access the vertex with global index "vertex" 
    //  * @param vertex: global index of the vertex
    // */
    // Vertex<d> vertex(int vertex) const {
    //     return vertices.at(vertex);
    // }

    // /**
    //  * @brief Access the element with global index "element" 
    //  * @param element: global index of the element
    // */
    // Quad<d> quad(int element) const {
    //     return quads.at(element);
    // }

    // Iterators for iterating over the quadrilaterals in the mesh
    using quad_iterator = typename std::list<Quad<d>>::const_iterator;  

    quad_iterator begin() const {
        return quads.begin();
    }

    quad_iterator end() const {
        return quads.end();
    }

};



struct Mesh1D : public Mesh<1> {

    Mesh1D(const double a, const double b, const int n);
};




#endif