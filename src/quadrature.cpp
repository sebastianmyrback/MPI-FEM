#include "../include/quadrature/quadrature.hpp"

// typedef quadrature::QuadratureRule<1>::QP QP;

// static std::vector<QP> qp_trapezoidal_1d = {
//     QP(Rd<1>(0.0), 0.5),
//     QP(Rd<1>(1.0), 0.5)
// };

// const quadrature::QuadratureRule<1> Trapezoidal1D(2, qp_trapezoidal_1d);

namespace quadrature {

    const QuadratureRule<1> midpoint(1, {
        QuadraturePoint<1>(Point<1>({0.5}), 1.)
        });

    const QuadratureRule<1> trapezoidal(2, {
        QuadraturePoint<1>(Point<1>({0.0}), 0.5), 
        QuadraturePoint<1>(Point<1>({1.0}), 0.5)
        });

    const QuadratureRule<1> simpson(3, {
        QuadraturePoint<1>(Point<1>({0.0}), 1./3), 
        QuadraturePoint<1>(Point<1>({0.5}), 1./3), 
        QuadraturePoint<1>(Point<1>({1.0}), 1./3)
        });

    const QuadratureRule<1> gauss_lobatto6(6, {
        QuadraturePoint<1>(Point<1>({0.0}), 0.03333333333333333), 
        QuadraturePoint<1>(Point<1>({0.11747233803526763}), 0.1892374781489235), 
        QuadraturePoint<1>(Point<1>({0.3573842417596774}), 0.2774291885177432),
        QuadraturePoint<1>(Point<1>({0.6426157582403226}), 0.2774291885177432), 
        QuadraturePoint<1>(Point<1>({0.8825276619647324}), 0.1892374781489235), 
        QuadraturePoint<1>(Point<1>({1.0}), 0.03333333333333333)
        });

} // namespace quadrature