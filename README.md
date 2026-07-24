# POLYNOMIAL-APPROXIMATION-USING-METHOD-OF-LEAST-SQUARES
POLYNOMIAL REGRESSION
# Polynomial Approximator using Householder QR Factorization

## Overview

This project implements polynomial approximation in C++ using the Least Squares Method and Householder QR Factorization. Given a set of data points stored in a text file, the program automatically determines the polynomial that best fits the data by constructing a Vandermonde matrix and solving the resulting least-squares problem. The polynomial degree is increased iteratively until the approximation satisfies the specified tolerance or the maximum permissible degree is reached.

The primary objective of this project is to demonstrate the implementation of fundamental numerical linear algebra algorithms from scratch without relying on external numerical libraries such as Eigen, LAPACK, or Armadillo. Every major component, including the matrix class, Householder QR factorization, back substitution, Vandermonde matrix generation, and residual computation, has been implemented manually in C++.

## Mathematical Background

Given a collection of data points ((x_i, y_i)), the program constructs a Vandermonde matrix corresponding to a polynomial of a chosen degree. The least-squares problem is formulated as (Ac \approx b), where (A) is the Vandermonde matrix, (c) is the vector of polynomial coefficients, and (b) contains the observed values. The matrix is factorized using Householder QR Factorization, after which the system (Rc = Q^Tb) is solved using back substitution to obtain the polynomial coefficients. The quality of the approximation is measured using the Euclidean norm of the residual vector.

## Project Structure

The project consists of a custom Matrix class (`MATRIX.cpp` and `matrix.h`), numerical utility functions (`NLA_FUNCTIONS.cpp` and `NLA_FUNCTIONS.h`), QR factorization routines (`factorizations.cpp` and `factorizations.h`), and an input dataset stored in `test.txt`. All source files are kept within the same project directory. The Visual Studio Code build configuration is provided through `tasks.json`, allowing the entire project to be compiled with a single build command.

## Input Format

The input dataset is read from `test.txt`. Each line of the file contains two space-separated values representing the x-coordinate and the corresponding y-coordinate of a data point. No header row is required. A typical input file is of the form:

```text
0 1
1 6
2 63
3 364
4 1365
5 3906
```

## Running the Project

Place all project files inside the same folder and update `test.txt` with the desired dataset before execution. Open the folder in Visual Studio Code and build the project using the configured build task (`Ctrl + Shift + B`). Once compilation is complete, execute the generated program. The application reads the dataset, performs polynomial approximation, and prints the best polynomial degree, residual norm, and the final polynomial equation.

## Features

The project provides automatic construction of the Vandermonde matrix, Householder QR factorization, back substitution for solving upper triangular systems, automatic polynomial degree selection based on the residual norm, and formatted printing of the resulting polynomial. The implementation is entirely based on the C++ Standard Library and does not depend on any external numerical computing packages.

## Output

After processing the dataset, the program reports the polynomial that best approximates the given data along with its degree and residual norm. A sample output is shown below:

```text
Best Polynomial Degree : 5

Residual Norm : 1.14216e-09

Best Fit Polynomial

P(x) = 1 + x + x² + x³ + x⁴ + x⁵
```

## Future Improvements

Future extensions may include support for additional basis functions such as Chebyshev and Legendre polynomials, implementation of weighted least squares, automatic model selection using statistical criteria, visualization of the fitted curve, and support for noisy datasets through regularization techniques.

## Author

Shyam Pithwa

Hackathon Project – Polynomial Approximation using Householder QR Factorization
