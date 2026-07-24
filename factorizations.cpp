#include "factorizations.h"
#include "NLA_FUNCTIONS.h"
#include <iostream>
#include <algorithm>

using namespace std;



Householder_QR::Householder_QR(const Matrix &A)
    : R(A), Q(identity(A.getRows()))
    {
        

        int n = A.getCols();
        int m = A.getRows();

        for (int i = 0; i < std::min(m - 1, n); i++)
        {  
            Matrix H = identity(m);
            Matrix A_sub = Sub_Matrix(R, i, m - 1, i, n - 1);


            vector<double> x = access_colunm(A_sub, 0);
            vector<double> e(x.size(), 0.0);
            e[0] = 1.0;

            double y = tnorm(x);
            double z = (x[0] >= 0) ? 1.0 : -1.0;

            vector<double> v = VS_Multiplication(e, y * z);
            v = Vec_Add(v, x);

            double v1 = tnorm(v);
            if (v1 == 0.0) continue;

            v = VS_Multiplication(v, 1.0 / v1);
            Matrix H_small=identity(v.size());
            H_small=H_small -(OuterProduct(v,v)*2.0) ;
            
            Matrix temp = VM_multiplication(v, A_sub);
            Matrix update = OuterProduct(v, temp) * 2.0;
            InsertSubMatrix(H,H_small,i,i);
            Q = Q*H;
            A_sub = A_sub - update;
            InsertSubMatrix(R, A_sub, i, i);
        }
    }
