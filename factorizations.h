#pragma once
#ifndef HOUSEHOLDER_QR_H
#define HOUSEHOLDER_QR_H

#include "matrix.h"

class Householder_QR
{
public:
    Matrix R;
    Matrix Q;
    Householder_QR(const Matrix &A);

    const Matrix& getQ() const
    {
        return Q;
    }

    const Matrix& getR() const
    {
        return R;
    }
};
vector<double> back_substitution(const Matrix &R, const Matrix &b);
vector<double> Polynomial_Approximator(const vector<double>& x,const vector<double> &y);
void printPolynomial(const vector<double>& coeff);
void readDataPoints(const char *fileName,vector<double> &x,vector<double> &y);
#endif