#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "mesh.hpp"

template<int d>
class QuadraturePoint {
public:

    const Rd<d> node;
    const double weight;

    QuadraturePoint(const Rd<d> &_node, const double _weight) : node(_node), weight(_weight) {}

};

template<int d>
class QuadratureRule {

public:
    
    typedef QuadraturePoint<d> QP;
    const int n;    // number of integratio points
    
private:
    
    std::vector<QP> points;

public:
    
    QuadratureRule(const int _n, std::vector<QP> qpoints) : n(_n), points(qpoints) {
        assert(n == qpoints.size());
    }
    
    const QP & operator[](int i) const {return points.at(i);}

};

typedef QuadratureRule<1> QuadratureRule1D;
extern const QuadratureRule1D Trapezoidal1D;


struct Quadrature1D {
    typedef QuadratureRule1D QuadratureRule;
    
    const QuadratureRule &set_quadrature_rule(const std::string &rule) {
        if (rule == "trapezoidal") {
            return Trapezoidal1D;
        }
        else return Trapezoidal1D;
    }   

};



template<int dim> struct TypeQuadrature {};
template<> struct TypeQuadrature<1> {
    typedef Quadrature1D Quadrature;
    typedef QuadraturePoint<1> QP;    
};


#endif // QUADRATURE_HPP