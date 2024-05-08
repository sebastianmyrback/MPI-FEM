#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "mesh.hpp"



template<int d>
class QuadraturePoint {
public:

    const Point<d> node;
    const double weight;

    QuadraturePoint(const Point<d> &_node, const double _weight) : node(_node), weight(_weight) {}

};

template<int d>
class QuadratureRule {

public:
    
    typedef QuadraturePoint<d> QP;      
    const int n;                        // number of quadrature points
    
private:
    
    std::vector<QP> points;

public:
    
    QuadratureRule(const int _n, std::vector<QP> qpoints) : n(_n), points(qpoints) {
        assert(n == qpoints.size());
    }
    
    const QP & operator[](int i) const {
        if (i >= n) {
            std::cerr << "QuadratureRule: index out of bounds" << std::endl;
            exit(1);
        }
        return points[i];
    }

};

#endif // QUADRATURE_HPP