// Point.tpp

template<int dim>
Point<dim> Point<dim>::operator+(const Point<dim> &other) const {
    Point<dim> result;
    for (std::size_t i = 0; i < dim; i++) {
        result.components[i] = this->components[i] + other.components[i];
    }
    return result;
}

template<int dim>
Point<dim> Point<dim>::operator-(const Point<dim> &other) const {
    Point<dim> result;
    for (std::size_t i = 0; i < dim; i++) {
        result.components[i] = this->components[i] - other.components[i];
    }
    return result;
}

template<int dim>
Point<dim> Point<dim>::operator*(double scalar) const {
    Point<dim> result;
    for (std::size_t i = 0; i < dim; i++) {
        result.components[i] = this->components[i] * scalar;
    }
    return result;
}

template<int dim>
Point<dim> Point<dim>::operator/(double scalar) const {
    assert(scalar != 0);
    Point<dim> result;
    for (std::size_t i = 0; i < dim; i++) {
        result.components[i] = this->components[i] / scalar;
    }
    return result;
}


template<int dim>
std::ostream &operator<<(std::ostream &os, const Point<dim> &p) {
    for (std::size_t i = 0; i < dim; i++) {
        os << p.components[i] << ' ';
    }
    return os;
}
