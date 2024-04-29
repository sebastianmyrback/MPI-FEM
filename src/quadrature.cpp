#include "quadrature.hpp"

typedef QuadratureRule<1>::QP QP;

static std::vector<QP> qp_trapezoidal_1d = {
    QP(Rd<1>(0.0), 0.5),
    QP(Rd<1>(1.0), 0.5)
};

const QuadratureRule<1> Trapezoidal1D(2, qp_trapezoidal_1d);
