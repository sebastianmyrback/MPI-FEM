#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <array>
#include <memory>
#include <iostream>

template<int d>
struct Rd {
    std::vector<double> x;

    // Default constructor
    Rd() : x(d, 0.0) {}

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
struct vertex : public Rd<d> {
    int vertex_label, idx;      // label and index

    // Default constructor
    vertex() : Rd<d>(), vertex_label(0) {}

    // Constructor that takes a double as an argument
    vertex(const double &s) : Rd<d>(std::vector<double>(d, s)), vertex_label(0) {}

    // Constructor that takes an Rd and an int as arguments
    vertex(const Rd<d> &P, int r = 0) : Rd<d>(P), vertex_label(r) {}

};

template<int d>
struct element {

    // Quadrilateral elements have 2^d nodes
    static const int n_vertices = 1 << d;

    // Pointer to the vertices of the element
    std::vector<std::shared_ptr<vertex<d>>> elem_vertices;

    const vertex<d> & operator()(int i) const {return *(elem_vertices.at(i));}

    double measure;

};


class mesh1d {
public :
    
    typedef vertex<1> vertex;
    typedef element<1> element;
    typedef Rd<1> Rd;

    static const int D = 1;

    int nv, nk;                     // number of nodes and elements
    double h;                       // typical mesh size
    std::vector<vertex> mesh_vertices;        // array of vertices
    std::vector<element> elements;        // array of elements

    mesh1d(double a, double b, int N);

    const vertex & operator()(int i) const {return mesh_vertices.at(i);} 
    const element & operator[](int i) const {return elements.at(i);} 

    //void build_mesh(double a, double b, int N);
  
};


class mesh2d {
public :
  typedef vertex<2> vertex;
  typedef element<2> element;

  static const int D = 2;
  
  int nv, nt;
  double area;
  std::vector<vertex> vertices;        // array of vertices
  std::vector<element> elements;        // array of elements
  
  //void BuildMesh(double a, double b, int N);
  
};



#endif