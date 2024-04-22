#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <fstream>


// 2D vector structure
class R2 {
public:
    double x, y;
    // Constructor
    R2(double _x, double _y) : x(_x), y(_y) {}
    R2(double s) : x(s), y(s) {}
    R2() : x(0), y(0) {}

    R2 &operator=(const R2 &other)
    {
        x = other.x;
        y = other.y;
        return *this;
    }

    R2 &operator=(const double &scalar)
    {
        x = scalar;
        y = scalar;
        return *this;
    }

    // Operator overloading
    double operator[](int i) const
    {
        return i == 0 ? x : y;
    }

    double operator()(int i) const
    {
        return i == 0 ? x : y;
    }

    R2 operator+(const R2 &other) const
    {
        return R2(x + other.x, y + other.y);
    }

    R2 operator-(const R2 &other) const
    {
        return R2(x - other.x, y - other.y);
    }

    R2 operator*(double scalar) const
    {
        return R2(x * scalar, y * scalar);
    }

    R2 operator/(double scalar) const
    {
        return R2(x / scalar, y / scalar);
    }

    R2 &operator+=(const R2 &other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    R2 &operator-=(const R2 &other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    R2 &operator*=(double scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    R2 &operator/=(double scalar)
    {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    bool operator==(const R2 &other) const
    {
        return x == other.x && y == other.y;
    }

    bool operator!=(const R2 &other) const
    {
        return x != other.x || y != other.y;
    }


    friend std::ostream &operator<<(std::ostream &f, const R2 &P) {
        f << P.x << ' ' << P.y;
        return f;
    }
    friend std::istream &operator>>(std::istream &f, R2 &P) {
        f >> P.x >> P.y;
        return f;
    }

};



class Vertex : public R2 {

public:

    int vertex_label;

    Vertex() : R2(), vertex_label(0) {};
    Vertex(const double &s) : R2(s), vertex_label(0) {};
    Vertex(const R2 &P, int r = 0) : R2(P), vertex_label(r) {};

};



#endif // VERTEX_HPP
