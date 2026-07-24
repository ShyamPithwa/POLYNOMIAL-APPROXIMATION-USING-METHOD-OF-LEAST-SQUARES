#pragma once
#ifndef __EE_242_Project_2__matrix__
#define __EE_242_Project_2__matrix__

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include<iomanip>
#include <cmath>

using std::tuple;
using std::vector;
class Matrix
{
private:
    unsigned m_rowSize;
    unsigned m_colSize;
    vector<vector<double>> m_matrix;

public:
    Matrix(unsigned, unsigned, double);
    Matrix(const char *);
    Matrix(const Matrix &);
    ~Matrix();

    // Matrix Operations
    Matrix operator+(Matrix &);
    Matrix operator-(const Matrix &) const;
    Matrix operator*(Matrix &);
    Matrix transpose();

    void input();

    // Scalar Operations
    Matrix operator+(double);
    Matrix operator-(double);
    Matrix operator*(double);
    Matrix operator/(double);

    // Aesthetic Methods
    double &operator()(const unsigned &rowNo, const unsigned &colNo);
    const double &operator()(const unsigned &rowNo, const unsigned &colNo) const;

    void print() const;
    Matrix clean(double eps = 1e-12) const;
    void prettyPrint(int precision = 6, double eps = 1e-12) const;

    unsigned getRows() const;
    unsigned getCols() const;

    // Power Iteration
    tuple<Matrix, double, int> powerIter(unsigned, double);

    // Deflation
    Matrix deflation(Matrix &, double &);
};
#endif