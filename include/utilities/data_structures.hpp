#pragma once 

#include <iostream>
#include <map>
#include <vector>

class DenseMatrix : public std::vector<std::vector<double>> 
{

    int n, m;   // rows, columns

public:
    
    DenseMatrix(int _n, int _m) : n(_n), m(_m) {
        this->assign(n, std::vector<double>(m, 0.0));
    }

    DenseMatrix(int _n, int _m, const double s) : n(_n), m(_m) {
        this->assign(n, std::vector<double>(m, s));
    }

    DenseMatrix(const DenseMatrix &A) : n(A.n), m(A.m), std::vector<std::vector<double>>(A) {}

    DenseMatrix(DenseMatrix &&A) : n(A.n), m(A.m), std::vector<std::vector<double>>(std::move(A)) {}

    DenseMatrix & operator=(const DenseMatrix &A) {
        n = A.n;
        m = A.m;
        std::vector<std::vector<double>>::operator=(A);
        return *this;
    }

    DenseMatrix & operator=(DenseMatrix &&A) {
        n = A.n;
        m = A.m;
        std::vector<std::vector<double>>::operator=(std::move(A));
        return *this;
    }

    double & operator()(int i, int j) {
        // if (i >= n || j >= m) {
        //     std::cerr << "DenseMatrix::operator(): index out of bounds" << std::endl;
        //     exit(1);
        // }
        return (*this)[i][j];
    }

    const double & operator()(int i, int j) const {
        // if (i >= n || j >= m) {
        //     std::cerr << "DenseMatrix::operator(): index out of bounds" << std::endl;
        //     exit(1);
        // }
        return (*this)[i][j];
    }

    
    void reinit(double val) {
        for (auto & row : (*this)) {
            std::fill(row.begin(), row.end(), val);
        }
    }

    int rows() const { return n; }

    int cols() const { return m; }

    void print() const {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                std::cout << (*this)(i, j) << ' ';
            }
            std::cout << std::endl;
        }
    }

};


class SparseMatrix : public std::map<std::pair<int, int>, double>
{

public:

    SparseMatrix() {}

    SparseMatrix(const SparseMatrix &A) : std::map<std::pair<int, int>, double>(A) {}

    SparseMatrix(SparseMatrix &&A) : std::map<std::pair<int, int>, double>(std::move(A)) {}

    SparseMatrix(const std::map<std::pair<int, int>, double> &A) : std::map<std::pair<int, int>, double>(A) {}

    SparseMatrix & operator=(const SparseMatrix &A) {
        (*this) = A;
        return *this;
    }

    SparseMatrix & operator=(SparseMatrix &&A) {
        (*this) = std::move(A);
        return *this;
    }

    void add(int i, int j, double val) {
        (*this)[std::make_pair(i, j)] += val;
    }

    void insert(int i, int j, double val) {
        (*this)[std::make_pair(i, j)] = val;
    }

    double operator()(int i, int j) const {
        auto it = (*this).find(std::make_pair(i, j));
        if (it != (*this).end()) {
            return it->second;
        }
        return 0.0;
    }

    void print() const {
     
        for (const auto & [indices, value] : (*this)) {
     
            int row = indices.first;
            int col = indices.second;

            std::cout << "(" << row << ", " << col << ") = " << value << std::endl;

        }
    }

};


class Vector : public std::vector<double> 
{

public:

    Vector() {}

    Vector(int n) {
        this->assign(n, 0.);
    }

    Vector(int n, const double s) {
        this->assign(n, s);
    }

    Vector(const Vector &v) : std::vector<double>(v) {} 

    Vector(Vector &&v) : std::vector<double>(std::move(v)) {}

    Vector(const std::vector<double> &v) : std::vector<double>(v) {}

    Vector & operator=(const Vector &v) {
        std::vector<double>::operator=(v);
        return *this;
    }

    Vector & operator=(Vector &&v) {   
        std::vector<double>::operator=(std::move(v));
        return *this;
    }

    double & operator()(int i) {
        if (i >= (*this).size()) {
            std::cerr << "Vector::operator(): index out of bounds" << std::endl;
            exit(1);
        }
        return (*this)[i];
    }

    const double & operator()(int i) const {
        if (i >= (*this).size()) {
            std::cerr << "Vector::operator(): index out of bounds" << std::endl;
            exit(1);
        }
        return (*this)[i];
    }

    // arithmetic operations
    Vector operator+(const Vector &v) const {
        if ((*this).size() != v.size()) {
            std::cerr << "Vector::operator+: dimension mismatch" << std::endl;
            exit(1);
        }
        Vector res(v.size());
        for (int i = 0; i < v.size(); i++) {
            res(i) = (*this)[i] + v(i);
        }
        return res;
    }

    Vector operator-(const Vector &v) const {
        if ((*this).size() != v.size()) {
            std::cerr << "Vector::operator-: dimension mismatch" << std::endl;
            exit(1);
        }
        Vector res(v.size());
        for (int i = 0; i < v.size(); i++) {
            res(i) = (*this)[i] - v(i);
        }
        return res;
    }

    Vector operator*(double scalar) const {
        Vector res((*this).size());
        for (int i = 0; i < (*this).size(); i++) {
            res(i) = (*this)[i] * scalar;
        }
        return res;
    }

    Vector operator/(double scalar) const {
        Vector res((*this).size());
        for (int i = 0; i < (*this).size(); i++) {
            res(i) = (*this)[i] / scalar;
        }
        return res;
    }

    Vector operator+=(const Vector &v) {
        if ((*this).size() != v.size()) {
            std::cerr << "Vector::operator+=: dimension mismatch" << std::endl;
            exit(1);
        }
        for (int i = 0; i < v.size(); i++) {
            (*this)[i] += v(i);
        }
        return *this;
    }

    Vector operator-=(const Vector &v) {
        if ((*this).size() != v.size()) {
            std::cerr << "Vector::operator-=: dimension mismatch" << std::endl;
            exit(1);
        }
        for (int i = 0; i < v.size(); i++) {
            (*this)[i] -= v(i);
        }
        return *this;
    }

    // void init(const size_t n, const double val)
    // {
    //     this->assign(n, val);
    // }

    void reinit(double val) {
        std::fill((*this).begin(), (*this).end(), val);
    }

    // dot product
    double dot(const Vector &v) const {
        double res = 0;
        for (int i = 0; i < (*this).size(); i++) {
            res += (*this)[i] * v(i);
        }
        return res;
    }

    // Matrix-Vector product
    Vector matvec(const DenseMatrix &A) const {
        
        if ((*this).size() != A.cols()) {
            std::cerr << "Vector::matvec: dimension mismatch" << std::endl;
            exit(1);
        }
        Vector res(A.cols());

        for (int i = 0; i < A.cols(); i++) {
            for (int j = 0; j < (*this).size(); j++) {
                res(i) += A(j, i) * (*this)[j];
            }
        }
        return res;
    }

    Vector matvec(const SparseMatrix &A) const 
    {
        // if ((*this).size() != A.size()) {
        //     std::cerr << "Vector::matvec: dimension mismatch" << std::endl;
        //     exit(1);
        // }

        Vector res((*this).size());
        for (auto & [indices, value] : A) {
            int row = indices.first;
            int col = indices.second;

            res.at(row) += value * (*this)[col];
        }

        return res;
    }

    friend Vector operator*(double scalar, const Vector &v) 
    {
        return v * scalar;
    }


    friend Vector operator*(const class DenseMatrix &A, const Vector &v) 
    {
        return v.matvec(A);
    }

    friend Vector operator*(const class SparseMatrix &A, const Vector &v) 
    {
        return v.matvec(A);
    }

    void print() const {

        for (int i = 0; i < (*this).size(); i++) {
            std::cout << (*this)(i) << std::endl;
        }
    }

    };

