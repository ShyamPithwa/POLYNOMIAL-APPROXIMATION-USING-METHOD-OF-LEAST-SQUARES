#include "matrix.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
// Helper functions

// 2-norm for matrix
double tnorm(const Matrix &A)
{
    if (A.getCols() != 1)
    {
        throw invalid_argument("tnorm() expects a column vector.");
    }

    double sum = 0.0;

    for (unsigned int i = 0; i < A.getRows(); i++)
    {
        sum += A(i, 0) * A(i, 0);
    }

    return sqrt(sum);
}

// 2-NORM for vector
double tnorm(const vector<double> &a)
{
  double n = a.size();
  double norm, s = 0.0;
  for (double i = 0; i < n; i++)
  {
    s += pow(a[i], 2);
  }
  norm = sqrt(s);
  return norm;
}

// 2 norm for Matrix colunm
double tnorm(const Matrix &A, int k)
{
  int n = A.getCols();
  double norm, s = 0.0;
  for (int i = 0; i < n; i++)
  {
    s += pow(A(i, k), 2);
  }
  norm = sqrt(s);
  return norm;
}

Matrix identity(int n)
{
  Matrix I(n, n, 0.0);

  for (int i = 0; i < n; i++)
    I(i, i) = 1.0;

  return I;
}

// Function to access colunm
vector<double> access_colunm(const Matrix &A, int k)
{
  int c = A.getRows();
  vector<double> colunm(c);
  for (int i = 0; i < c; i++)
  {
    colunm[i] = A(i, k);
  }
  return colunm;
}

// vector matrix outerproduct
Matrix OuterProduct(const vector<double> &v, const Matrix &row)
{
  // row must be a 1 × n matrix
  if (row.getRows() != 1)
  {
    throw invalid_argument("Second argument must be a row vector.");
  }

  Matrix result(v.size(), row.getCols(), 0.0);

  for (unsigned i = 0; i < v.size(); i++)
  {
    for (unsigned j = 0; j < row.getCols(); j++)
    {
      result(i, j) = v[i] * row(0, j);
    }
  }

  return result;
}

// vector vector outer product
Matrix OuterProduct(const vector<double> &a, const vector<double> &b)
{
  Matrix result(a.size(), b.size(), 0.0);

  for (unsigned i = 0; i < a.size(); i++)
  {
    for (unsigned j = 0; j < b.size(); j++)
    {
      result(i, j) = a[i] * b[j];
    }
  }

  return result;
}

// For Multiplication of vector with matrix
Matrix VM_multiplication(const vector<double> &v, const Matrix &A)
{
  if (v.size() != A.getRows())
    throw invalid_argument("Dimension mismatch.");

  Matrix result(1, A.getCols(), 0.0);
  for (int j = 0; j < A.getCols(); j++)
  {
    double sum = 0.0;

    for (int i = 0; i < A.getRows(); i++)
    {
      sum += v[i] * A(i, j);
    }

    result(0, j) = sum;
  }

  return result;
}

// for multiplication of matrix with vector
Matrix MV_multiplication(const Matrix &A, const vector<double> &v)
{
  if (v.size() != A.getCols())
    throw invalid_argument("Dimension mismatch.");

  Matrix result(A.getRows(), 1, 0.0);

  for (int i = 0; i < A.getRows(); i++)
  {
    double sum = 0.0;

    for (int j = 0; j < A.getCols(); j++)
    {
      sum += A(i, j) * v[j];
    }

    result(i, 0) = sum;
  }

  return result;
}
// For Multiplication of scalar to whole vector
vector<double> VS_Multiplication(const vector<double> &v, double k)
{
  int n = v.size();
  vector<double> result(n); // Create a clean new vector
  for (int i = 0; i < n; i++)
  {
    result[i] = v[i] * k;
  }
  return result; // Return the new vector, leaving 'v' untouched!
}

// Dot/inner product
double Dot_product(const vector<double> Q, const vector<double> P)
{

  if (Q.size() == P.size())
  {
    int n = Q.size();
    double r = 0;
    for (int i = 0; i < n; i++)
    {
      r += P[i] * Q[i];
    }
    return r;
  }

  else
  {
    cout << "INVALID OPERATION" << endl;
    return -1;
  }
}

// Vector Subtraction
vector<double> Vec_Sub(const vector<double> A, const vector<double> B)
{
  if (A.size() == B.size())
  {
    int n = A.size();
    vector<double> C(n);
    for (int i = 0; i < n; i++)
    {
      C[i] = A[i] - B[i];
    }
    return C;
  }
  else
    return {};
}

// Vector addition
vector<double> Vec_Add(const vector<double> &A, const vector<double> &B)
{
  if (A.size() == B.size())
  {
    vector<double> C(A.size());
    for (int i = 0; i < A.size(); i++)
    {
      C[i] = A[i] + B[i];
    }
    return C;
  }

  else
  {
    throw invalid_argument("VECTORS ARE OF DIFRENT DIMENSIONS!!!");
  }
}

// To update colunms
void set_column(Matrix &A, int k, const vector<double> &col)
{
  if (col.size() != A.getRows())
  {
    cout << "Invalid column size" << endl;
    return;
  }

  for (int i = 0; i < A.getRows(); i++)
  {
    A(i, k) = col[i];
  }
}

// Modified Gram schmidt
Matrix MG_Schmidt(const Matrix &A)
{
  int m = A.getCols();
  int n = A.getRows();
  Matrix V = A;
  for (int i = 0; i < m; i++)
  {
    vector<double> v = access_colunm(V, i);
    double r = tnorm(v);
    if (r != 0)
    {
      vector<double> q = VS_Multiplication(v, 1.0 / r);
      set_column(V, i, q);
      for (int j = i + 1; j < m; j++)
      {
        vector<double> vj = access_colunm(V, j);
        double z = Dot_product(q, vj);
        vector<double> proj = VS_Multiplication(q, z);
        vj = Vec_Sub(vj, proj);
        set_column(V, j, vj);
      }
    }
    else
    {
      cout << "NOT POSSIBLE" << endl;
    }
  }
  return V;
}

// Sub matrix function
Matrix Sub_Matrix(const Matrix &A, int i1, int i2, int j3, int j4)
{
  if (i1 < 0 || i2 >= static_cast<int>(A.getRows()) ||
      j3 < 0 || j4 >= static_cast<int>(A.getCols()) ||
      i1 > i2 || j3 > j4)
  {
    throw invalid_argument("ENTER CORRECT VALUES TO OBTAIN SUBMATRIX!!!");
  }

  Matrix c(i2 - i1 + 1, j4 - j3 + 1, 0.0);
  for (int i = i1, r = 0; i <= i2; i++, r++)
  {
    for (int j = j3, ccol = 0; j <= j4; j++, ccol++)
    {
      c(r, ccol) = A(i, j);
    }
  }
  return c;
}

void InsertSubMatrix(Matrix &A, const Matrix &sub, unsigned startRow, unsigned startCol)
{
  if (startRow + sub.getRows() > A.getRows() ||
      startCol + sub.getCols() > A.getCols())
  {
    throw invalid_argument("Submatrix exceeds matrix dimensions.");
  }

  for (unsigned i = 0; i < sub.getRows(); i++)
  {
    for (unsigned j = 0; j < sub.getCols(); j++)
    {
      A(startRow + i, startCol + j) = sub(i, j);
    }
  }
}

// To construct vandermonde matrix
// Add degree by considering highest degree of polynomial
Matrix vandermonde(const vector<double> &x, int degree)
{
  Matrix V(x.size(), degree + 1, 0.0);

  for (int i = 0; i < x.size(); i++)
  {
    double value = 1.0;
    for (int j = 0; j <= degree; j++)
    {
      V(i, j) = value;
      value *= x[i];
    }
  }

  return V;
}

// vector to matrix
Matrix vectorToMatrix(const vector<double> &v)
{
  Matrix M(v.size(), 1, 0.0);

  for (int i = 0; i < v.size(); i++)
  {
    M(i, 0) = v[i];
  }

  return M;
}