#include "eigne_approximators.h"

double RQI(const Matrix &A)
{
    int m = A.getRows();
    int n = A.getCols();
    int iterations = 0;
    if (m != n)
    {
        throw std::invalid_argument("Matrix must be square for RQI.");
    }
    Matrix v(m, 1, 1);
    v = v * (1.0 / tnorm(access_colunm(v, 0)));
    Matrix p = (v.transpose()) * (A) * (v);
    double pe = p(0, 0);
    double ce = 0.0;
    Matrix I = identity(m);
    while (true)
    {
        if(iterations==m)
        {
            std::cout << "Maximum iterations reached. The method may not have converged." << std::endl;
            break;
        }
        Matrix B = A - (I * pe);
        Householder_QR qr(B);
        Matrix R1 = qr.R;
        Matrix Qt = qr.Q.transpose();
        Matrix b = Qt * v;
        Matrix rhs = Sub_Matrix(b, 0, m - 1, 0, 0);
        v = vectorToMatrix(back_substitution(R1, rhs));
        v = v * (1.0 / tnorm(access_colunm(v, 0)));
        Matrix temp = (v.transpose()) * (A) * (v);
        ce = temp(0, 0);
        Matrix res = (A * v) - (v * ce);
        double res_norm = tnorm(access_colunm(res, 0));
        if (res_norm < 1e-12){
            break;
        }
        else
        {
            pe = ce;
        }
        iterations++;
    }
    std::cout << "Converged in " << iterations << " iterations." << std::endl;
    return ce;
}
