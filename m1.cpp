#include <iostream>
#include "factorizations.h"
#include "matrix.h"
#include "NLA_FUNCTIONS.h"
using namespace std;

int main()
{
  cout<<"\nPOLYNOMIAL APPROXIMATOR\n" ;
  cout<<"\n ENTER DATA SET(txt file)\n";
  vector<double>x,y;
  readDataPoints("test.txt",x,y);
  vector<double> c= Polynomial_Approximator(x,y);
  printPolynomial(c);
}
