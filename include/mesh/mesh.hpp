#pragma once

#include <vector>
#include <array>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
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

    size_t size() const noexcept { return dim; }

    double& operator[](size_t i) {
        assert(i < dim);
        return components[i];
    }

    const double& operator[](size_t i) const {
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
        for (size_t i = 0; i < dim; i++) {
            os << p.components[i] << ' ';
        }
        return os;
    }

};



template<int dim>
class Vertex : public Point<dim> {

        
    size_t vertex_label, glb_idx;      // label and index

public:
    inline size_t boundary_label() const noexcept {return vertex_label;}
    inline size_t global_index() const noexcept {return glb_idx;}

    Vertex() : Point<dim>() {assert(false);}

    Vertex(const Point<dim> &p) : Point<dim>(p), vertex_label(0), glb_idx(0) {}

    Vertex(const Point<dim> &p, const size_t idx, const size_t label) : Point<dim>(p), glb_idx(idx), vertex_label(label) {}
};


template<int dim>
struct Cell {

    static constexpr size_t n_verts_per_cell = 1 << dim;   // number of vertices in a quadrilateral cell = 2^d

    // Access the vertex with local index vertex 
    inline const Vertex<dim> &vertex(size_t vertex) const {
    
        size_t vertex_index = msh->get_cell_to_vertex()[index*n_verts_per_cell + vertex];

        return msh->get_vertices()[vertex_index];
    }

    // Map a point from the reference element xref to the physical element X
    inline void map_to_physical(const Point<dim> &xref, Point<dim> &x) const {
        
        x = vertex(0) + xref[0] * (vertex(1) - vertex(0));
    };

    Cell() : msh(nullptr), index(0), subdomain(0), measure(0.0) {assert(false);}
    Cell(Mesh<dim> *_msh, size_t i, double m) : msh(_msh), index(i), subdomain(0), measure(m) {}
    Cell(Mesh<dim> *_msh, size_t i, size_t d, double m) : msh(_msh), index(i), subdomain(d), measure(m) {}

    inline size_t get_index() const noexcept {return index;}
    inline size_t get_subdomain() const noexcept {return subdomain;}
    inline double get_measure() const noexcept {return measure;}
    
    void set_subdomain(size_t d) {subdomain = d;}

protected:

    const Mesh<dim> *msh;

    size_t index, subdomain;  // index of the cell in the mesh, subdomain (for parallelization)
    double measure;                // measure of the cell (length, area, volume, etc.)

};


template<int d>
struct Mesh {

protected:

    size_t nverts, ncells, nbe, nsubdomains;    // number of vertices, quadrilaterals and border elements in the mesh
    double h;                                   // typical mesh size

    std::vector<Cell<d>> cells;                 // list of cells in the mesh
    std::vector<Vertex<d>> vertices;            // list of vertices in the mesh
    std::vector<size_t> cell_to_vertex;         // map index from cell to vertex

public:

    constexpr static size_t dim = d;            // dimension of the mesh

    virtual void refine_mesh() = 0; 
    virtual void partition(const size_t n_mpi_processes) = 0;               // partition the mesh for parallelization   

    const size_t &get_nverts() const noexcept {return nverts;}
    const size_t &get_ncells() const noexcept {return ncells;}
    const size_t &get_nbe() const noexcept {return nbe;}
    const size_t &get_nsubdomains() const noexcept {return nsubdomains;}
    const double &get_h() const noexcept {return h;}

    virtual void mesh_info(const bool detailed = false) = 0;

    const std::vector<Vertex<d>> &get_vertices() const noexcept {return vertices;}
    const std::vector<size_t> &get_cell_to_vertex() const noexcept {return cell_to_vertex;}


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
    void refine_mesh() override;
    void partition(const size_t n_mpi_processes) override;
    void distribute_dofs();
    
    // const std::vector<std::vector<size_t>> &get_distribution() override;  
    // const std::map<int, std::set<int>> &get_shared_dofs() override;

    const std::vector<std::vector<size_t>> &get_distribution() {return dof_distribution;}
    const std::map<int, std::set<int>> &get_shared_dofs() {return shared_dofs;}

    void mesh_info(const bool detailed = false) override;

private:
    std::vector<std::vector<size_t>> dof_distribution;
    std::map<int, std::set<int>> shared_dofs;
    
};



#include "mesh.tpp"

