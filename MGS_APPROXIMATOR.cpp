#include<iostream>
#include "factorizations.h"
#include "NLA_FUNCTIONS.h"
#include <vector>
using namespace std;

int main(){
    vector<double> x ,y;
    readDataPoints("test.txt",x,y);
    printPolynomial(Polynomial_Approximator_MGS(x,y));
}