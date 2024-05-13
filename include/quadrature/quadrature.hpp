#pragma once

#include "../mesh/mesh.hpp"

namespace quadrature 
{

    template<int d>
    class QuadraturePoint 
    {
    public:

        const Point<d> node;
        const double weight;

        QuadraturePoint(const Point<d> &_node, const double _weight) : node(_node), weight(_weight) {}

    };

    template<int d>
    class QuadratureRule 
    {

    public:
        
        typedef QuadraturePoint<d> QP;      
        const int n;                        // number of quadrature points
        
    private:
        
        std::vector<QP> points;

    public:
        
        QuadratureRule(const int _n, std::vector<QP> qpoints) : n(_n), points(qpoints) {
            assert(n == qpoints.size());
        }

        QuadratureRule(const QuadratureRule &qr) : n(qr.n), points(qr.points) {}
        
        const QP & operator[](int i) const {
            if (i >= n) {
                std::cerr << "quadrature::QuadratureRule: index out of bounds" << std::endl;
                exit(1);
            }
            return points[i];
        }

    };


    extern const QuadratureRule<1> midpoint;
    extern const QuadratureRule<1> trapezoidal;
    extern const QuadratureRule<1> simpson;
    extern const QuadratureRule<1> gauss_lobatto6;


} // namespace quadrature