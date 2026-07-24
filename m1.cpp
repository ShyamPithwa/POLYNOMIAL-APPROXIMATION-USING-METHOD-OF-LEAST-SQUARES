#include <iostream>
#include "factorizations.h"
#include "matrix.h"
#include "NLA_FUNCTIONS.h"
using namespace std;

int main()
{
    int rows, cols;

    cout << "Householder QR Factorization\n";
    cout << "----------------------------\n";

    Matrix A("test.txt");

    Householder_QR qr(A);
    cout << "\nInput matrix A:\n";
    cout << "\nUpper-triangular matrix R:\n";
    qr.R.prettyPrint();
    qr.R.print();
    cout<<"\nOrthogonal Matrix Q: \n";
    qr.Q.print();
    qr.Q.prettyPrint();

    Matrix Z = (qr.Q)*(qr.R);
    cout<<"QR product matrix: "<<endl;
    Z.print();
    cout<<"Orignial matrix : "<<endl;
    A.print();

    cout<<"Orthogonality check for Q: "<<endl;
    Matrix k= qr.Q;
    Matrix p=qr.Q.transpose();
    Matrix i = k*p;
    i.print();
    i.prettyPrint();


    cout<<"Approximating/making both look clean:"<<endl;
    Z.prettyPrint();
    A.prettyPrint();
    return 0;
}
