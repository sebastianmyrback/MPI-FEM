#ifndef GNUPLOT_HPP_
#define GNUPLOT_HPP_

#include <string>
#include <fstream>
#include "mesh.hpp"

namespace gnuplot {

void save(const Mesh &Th, std::string filename = "Th.dat") {

    std::ofstream plot;
    plot.open(filename.c_str(), std::ofstream::out);
    const int nve = Th.get_element(0).n_vertices;
    for (int k = 0; k < Th.get_n_elements(); ++k) {
        for (int i = 0; i < nve; ++i) {
            plot << Th.get_element(k).get_vertex(i) << std::endl;
        }
        plot << Th.get_element(k).get_vertex(0) << std::endl;
        plot << std::endl;
        plot << std::endl;
    }
    plot.close();
}


} // namespace gnuplot

#endif
