#pragma once

#include <string>
#include <fstream>
#include <iomanip>
#include "../mesh/mesh.hpp"

namespace gnuplot {

template<int dim>
void write_cells(const Mesh<dim> &Th, std::string filename = "Th.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    
    for (auto cell = Th.cell_begin(); cell != Th.cell_end(); ++cell) {
        for (int i = 0; i < cell->n_verts_per_cell; ++i) {
            plot << cell->vertex(i) << std::endl;
        }
    }
    
    plot.close();
}

// void save(const mesh1d &Th, std::string filename = "Th.dat") {

//     std::ofstream plot;
//     plot.open(filename.c_str(), std::ofstream::out);
//     const int nve = Th.vertices[0].size();
//     for (int k = 0; k < Th.nk; ++k) {
//         for (int i = 0; i < nve; ++i) {
//             plot << *(Th.elements.at(k).vertices.at(i)) << std::endl;
//         }
//         plot << *(Th.elements.at(k).vertices.at(0)) << std::endl;
//         plot << std::endl;
//         plot << std::endl;
//     }
//     plot.close();
// }


// void save(const mesh2d &Th, std::string filename = "Th.dat") {

//     std::ofstream plot;
//     plot.open(filename.c_str(), std::ofstream::out);
//     const int nve = Th.get_element(0).n_vertices;
//     for (int k = 0; k < Th.get_n_elements(); ++k) {
//         for (int i = 0; i < nve; ++i) {
//             plot << Th.get_element(k).get_vertex(i) << std::endl;
//         }
//         plot << Th.get_element(k).get_vertex(0) << std::endl;
//         plot << std::endl;
//         plot << std::endl;
//     }
//     plot.close();
// }


} // namespace gnuplot


namespace matlab {
    
    void save(std::map<std::pair<int, int>, double> &dF, std::string filename) {
    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    plot << std::setprecision(16);

    for (auto q = dF.begin(); q != dF.end(); ++q) {
        plot << q->first.first + 1 << "\t" << q->first.second + 1 << "\t"
            << q->second << std::endl;
    }

    // for (std::map<std::pair<int, int>, double>::const_iterator q = dF.begin();
    //         q != dF.end(); ++q) {
    //     plot << q->first.first + 1 << "\t" << q->first.second + 1 << "\t"
    //         << q->second << std::endl;
    // }
    plot.close();
    }

    void save(std::vector<double> &rhs, std::string filename) {
        std::ofstream plot;
        plot.open(filename.c_str(), std::ofstream::out);
        plot << std::setprecision(16);
        for (int i = 0; i < rhs.size(); ++i) {
            plot << "\t" << rhs[i] << std::endl;
        }
        plot.close();
    }


};   // namespace matlab


