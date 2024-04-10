#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <array>
#include <fstream>


// 2D vector structure
class R2
{
public:
    double x, y;
    // Constructor
    R2(double _x, double _y) : x(_x), y(_y) {}
    R2() {}

    R2 &operator=(const R2 &other)
    {
        x = other.x;
        y = other.y;
        return *this;
    }

    R2 &operator=(double scalar)
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

    double Dot(const R2 &other) const
    {
        return x * other.x + y * other.y;
    }

    double Cross(const R2 &other) const
    {
        return x * other.y - y * other.x;
    }


    friend std::ostream &operator<<(std::ostream &os, const R2 &v)
    {
        os << "(" << v.x << ", " << v.y << ")";
        return os;
    }
};

// Quadrilateral element structure
class Element {

public:

    typedef R2 Vertex;
    typedef std::array<R2, 2> Edge;
    
    const int n_vertices = 4;
    const int n_edges = 4;

    // Constructor
    Element(std::array<Vertex, 4> _vertices, std::array<Edge, 4> _edges)
    {
        for (int i = 0; i < 4; i++)
        {
            vertices[i] = _vertices[i];
            edges[i] = _edges[i];
        }
    }

    Element()
    {
        for (int i = 0; i < 4; i++)
        {
            vertices[i] = 0.;
            edges[i][0] = 0.;
            edges[i][1] = 0.;
        }
    }

    // // Get vertex
    // const Vertex &GetVertex(unsigned int i) const
    // {
    //     return vertices[i];
    // }

    // // Get face
    // const Edge &GetEdge(unsigned int i) const
    // {
    //     return edges[i];
    // }

    // // Set vertex
    // void SetVertex(unsigned int i, const Vertex &vertex)
    // {
    //     vertices[i] = vertex;
    // }

    // // Set face
    // void SetFace(unsigned int i, const Edge &edge)
    // {
    //     edges[i] = edge;
    // }

    // Get area
    double get_area() const
    {
        return vertices[1].x * vertices[2].y; // base*height
    }


private:
    std::array<Vertex, 4> vertices;
    std::array<Edge, 4> edges;

};

#endif // VERTEX_HPP
