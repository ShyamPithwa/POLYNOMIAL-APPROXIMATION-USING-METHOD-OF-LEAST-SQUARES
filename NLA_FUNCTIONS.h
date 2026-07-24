#pragma once
#ifndef GRAMSCHMIDT_H
#define GRAMSCHMIDT_H

#include "matrix.h"

Matrix MG_Schmidt(Matrix &A);

Matrix identity(int n);
Matrix Sub_Matrix(const Matrix &A, int i1, int i2, int j3, int j4);
void InsertSubMatrix(Matrix& A,const Matrix& sub,unsigned startRow,unsigned startCol);

double tnorm(const vector<double>& a);
double tnorm(Matrix &A, int k);

vector<double> access_colunm(const Matrix &A, int k);
void set_column(const Matrix &A, int k, const vector<double>& col);

vector<double> VS_Multiplication(const vector<double>& v, double k);
double Dot_product(const vector<double>& Q,const vector<double>& P);
Matrix VM_multiplication(const vector<double>& v,const Matrix& A);
Matrix OuterProduct(const vector<double>& v, const Matrix& row);
Matrix OuterProduct(const vector<double> &a,const vector<double> &b);

vector<double> Vec_Sub(const vector<double>& A,const vector<double>& B);
vector<double> Vec_Add(const vector<double>& A,const vector<double>& B);

#endif