#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <assert.h>

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
struct vertex {
    
    Rd<d> x;                      // coordinates
    
    int vertex_label, glb_idx;      // label and index

};

template<int d>
struct base_element {

    // Quadrilateral elements have 2^d nodes
    static const int n_vertices = 1 << d;

    // Pointers to the vertices of the element
    std::vector<std::shared_ptr<vertex<d>>> elem_vertices;

    const vertex<d> & operator()(int i) const {return *(elem_vertices.at(i));}

    double measure;

    virtual void operator()(const Rd<d> &xref, Rd<d> &x) const = 0;

};

// Template specialization for element

template<int d>
struct element : public base_element<d> {
    element() {assert(false);}
};

template<>
struct element<1> : public base_element<1> {
    element() {
        elem_vertices.resize(2);
    }

    void operator()(const Rd<1> &xref, Rd<1> &x) const {
        //x = elem_vertices[0]->x * (1.0 - xref[0]) + elem_vertices[1]->x * xref[0];
        

        x = elem_vertices[0]->x + xref[0] * (elem_vertices[1]->x - elem_vertices[0]->x);
        //Rd<1> x = elem_vertices[0]->x * (1.0 - xref[0]) + elem_vertices[1]->x * xref[0];
    }
};

template<int d>
class mesh {

public:
    
        typedef vertex<d> vert;
        typedef element<d> elem;
        typedef Rd<d> Rn;
        static const int D = d;

        int nv, nk;                             // number of nodes and elements
        double h;                               // typical mesh size
        std::vector<vert> mesh_vertices;        // array of vertices
        std::vector<elem> elements;             // array of elements
        std::vector<int> border_dofs;           // list of global border dofs

        mesh() : nv(0), nk(0) {}
    
        const vert & operator()(int i) const {return mesh_vertices.at(i);} 
        const elem & operator[](int i) const {return elements.at(i);} 
    
        //void build_mesh(double a, double b, int N);

};


class mesh1d : public mesh<1> {
public :
    
    mesh1d(double a, double b, int N);
  
};


class mesh2d {
public :
  
    mesh2d(double x0, double y0, double xend, double yend, int nx, int ny);

};



#endif