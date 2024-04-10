// Implement the mesh class

#include <assert.h>
#include "mesh.hpp"
#include "element.hpp"

// // Constructor
// Mesh::Mesh(){}

// // Constructor
// Mesh::Mesh(int _nx, int _ny, double _x0, double _y0, double _lx, double _ly) : nx(_nx), ny(_ny), x0(_x0), y0(_y0), lx(_lx), ly(_ly) {

//     assert(n_vertices == 0 && n_elements == 0 && n_be == 0);

//     std::array<int, 4> indQ = {0, 1, 3, 2};

//     n_vertices = nx * ny;
//     n_elements = (nx - 1) * (ny - 1);
//     n_be = 2 * ((nx - 1) + (ny - 1));
//     const double hx = lx / (nx - 1);
//     const double hy = ly / (ny - 1);




//     //KN<int> iv(4), indT(4);

//     int jt = 0;
//     for (int j = 0; j < ny - 1; j++) {
//         for (int i = 0; i < nx - 1; i++) {

//             int id = 0;
//             for (int jj = j; jj < j + 2; ++jj) {
//                 for (int ii = i; ii < i + 2; ++ii) {

//                     int ivl  = ii + jj * nx; // index
//                     iv(id++) = ivl;

//                     vertices[ivl].x = ii * hx + orx;
//                     vertices[ivl].y = jj * hy + ory;
//                 }
//             }
//             for (int e = 0; e < 4; ++e) {
//                 indT(e) = iv(indQ[e]);
//             }
//             elements[jt++].set(vertices, indT, 0);
//         }
//     }

//     // create the for borders
//     int lab, k = 0;
//     for (int i = 0; i < nx - 1; ++i) {
//         indT(0) = i;
//         indT(1) = i + 1;
//         lab     = 1;
//         for (int j = 0; j < 2; ++j)
//             vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
//         borderelements[k++].set(vertices, indT, lab);
//     }
//     for (int i = 0; i < ny - 1; ++i) {
//         indT(0) = (i + 1) * nx - 1;
//         indT(1) = indT(0) + nx;
//         lab     = 2;
//         for (int j = 0; j < 2; ++j)
//             vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
//         borderelements[k++].set(vertices, indT, lab);
//     }
//     for (int i = 0; i < nx - 1; ++i) {
//         indT(0) = i + nx * (ny - 1);
//         indT(1) = indT(0) + 1;
//         lab     = 3;
//         for (int j = 0; j < 2; ++j)
//             vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
//         borderelements[k++].set(vertices, indT, lab);
//     }
//     for (int i = 0; i < ny - 1; ++i) {
//         indT(0) = i * nx;
//         indT(1) = indT(0) + nx;
//         lab     = 4;
//         for (int j = 0; j < 2; ++j)
//             vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
//         borderelements[k++].set(vertices, indT, lab);
//     }
    
// }