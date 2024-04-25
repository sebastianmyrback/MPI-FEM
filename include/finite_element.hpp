#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP


template<typename mesh>
class finite_element {

public:
    typedef typename mesh::element element;
    typedef typename mesh::vertex vertex;
    typedef typename mesh::Rd Rd;

    // Constructor
    finite_element() {}

    // Destructor
    ~finite_element() {}

    // Evaluate the basis functions
    virtual void eval_basis(const element &elem, const Rd &x, std::vector<double> &phi) const = 0;

    // Evaluate the gradient of the basis functions
    virtual void eval_basis_grad(const element &elem, const Rd &x, std::vector<Rd> &grad) const = 0;

};

#endif