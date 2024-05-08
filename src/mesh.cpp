#include "mesh.hpp"

// // template specialization for 1D
// template<>
// void Mesh<1>::build_mesh(const double a, const double b, const int n) {
//     vertices.clear();
//     quads.clear();
//     quad_to_vertex.clear();
//     border_dofs.clear();

//     // create vertices
//     for (int i = 0; i < n + 1; i++) {
//         vertices.push_back(Vertex<1>({a + i * (b - a) / n}));
//     }

//     // create quads
//     for (int i = 0; i < n; i++) {
//         quads.push_back(Quad<1>(i, vertices[i], vertices[i + 1]));
//     }

//     // create quad_to_vertex map
//     quad_to_vertex.resize(n);
//     for (int i = 0; i < n; i++) {
//         quad_to_vertex[i] = i;
//     }

//     // create border dofs
//     border_dofs.push_back(0);
//     border_dofs.push_back(n);
// }


Mesh1D::Mesh1D(
    const double a, 
    const double b, 
    const int n) 
    {

    nverts = n + 1;
    nquads = n;
    nbe = 2;

    h = (b - a) / nquads;

    vertices.resize(nverts);
    quad_to_vertex.resize(n_verts_per_quad * nquads);
    border_dofs.reserve(nbe);

    // Create inner vertices
    for (int i = 1; i < nverts - 1; i++) {
        const Rd<1> x(a + i*h);  
        vertices[i] = Vertex<1>(x, i, 0);
        
    }

    // Mark the boundary vertices
    vertices[0] = Vertex<1>(a, 0, 1);
    vertices[nverts-1] = Vertex<1>(b, nverts-1, 2);

    border_dofs.push_back(0);
    border_dofs.push_back(nverts-1);

    // Create elements
    int l = 0;
    for (int i = 0; i < nquads; i++) {

        for (int j = 0; j < n_verts_per_quad; j++) {
            quad_to_vertex[i * n_verts_per_quad + j] = i + j;
        }

        quads.push_back(Quad<1>(std::make_shared<Mesh1D>(*this), i, h));      
    
    }


    return;

}
