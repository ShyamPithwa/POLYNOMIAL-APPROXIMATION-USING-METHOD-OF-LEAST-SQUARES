#include "matrix.h"
using namespace std;

// Constructor for Any Matrix
Matrix::Matrix(unsigned rowSize, unsigned colSize, double initial)
{
    m_rowSize = rowSize;
    m_colSize = colSize;
    m_matrix.resize(rowSize);
    for (unsigned i = 0; i < m_matrix.size(); i++)
    {
        m_matrix[i].resize(colSize, initial);
    }
}

void Matrix::input()
{
    std::cout << "Enter matrix elements row-wise:\n";

    for (unsigned i = 0; i < m_rowSize; i++)
    {
        for (unsigned j = 0; j < m_colSize; j++)
        {
            std::cin >> m_matrix[i][j];
        }
    }
}

// Constructor for Given Matrix
Matrix::Matrix(const char *fileName)
{
    ifstream file_A(fileName);

    if (!file_A.is_open())
    {
        throw runtime_error("Failed to open file");
    }

    string line;
    vector<double> row;

    while (getline(file_A, line))
    {
        stringstream ss(line);

        row.clear();

        double element;

        while (ss >> element)
        {
            row.push_back(element);
        }

        if (!row.empty())
        {
            m_matrix.push_back(row);
        }
    }

    m_rowSize = m_matrix.size();

    if (m_rowSize > 0)
        m_colSize = m_matrix[0].size();
    else
        m_colSize = 0;

    file_A.close();
}

// Copy Constructor
Matrix::Matrix(const Matrix &B)
{
    this->m_colSize = B.getCols();
    this->m_rowSize = B.getRows();
    this->m_matrix = B.m_matrix;
}

Matrix::~Matrix()
{
}

// Addition of Two Matrices
Matrix Matrix::operator+(Matrix &B)
{
    Matrix sum(m_colSize, m_rowSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++)
    {
        for (j = 0; j < m_colSize; j++)
        {
            sum(i, j) = this->m_matrix[i][j] + B(i, j);
        }
    }
    return sum;
}

// Subtraction of Two Matrices
Matrix Matrix::operator-(const Matrix &B) const
{
    Matrix diff(m_rowSize, m_colSize, 0.0);

    for (unsigned i = 0; i < m_rowSize; i++)
    {
        for (unsigned j = 0; j < m_colSize; j++)
        {
            diff(i, j) = (*this)(i, j) - B(i, j);
        }
    }

    return diff;
}

// Multiplication of Two Matrices
Matrix Matrix::operator*(Matrix &B)
{
    Matrix multip(m_rowSize, B.getCols(), 0.0);
    if (m_colSize == B.getRows())
    {
        unsigned i, j, k;
        double temp = 0.0;
        for (i = 0; i < m_rowSize; i++)
        {
            for (j = 0; j < B.getCols(); j++)
            {
                temp = 0.0;
                for (k = 0; k < m_colSize; k++)
                {
                    temp += m_matrix[i][k] * B(k, j);
                }
                multip(i, j) = temp;
                // cout << multip(i,j) << " ";
            }
            // cout << endl;
        }
        return multip;
    }
    else
    {
        throw invalid_argument("Dimension mismatch");
    }
}

Matrix Matrix::operator*(const Matrix &B) const
{
    if (m_colSize != B.getRows())
    {
        throw invalid_argument("Dimension mismatch");
    }

    Matrix multip(m_rowSize, B.getCols(), 0.0);

    for (unsigned i = 0; i < m_rowSize; i++)
    {
        for (unsigned j = 0; j < B.getCols(); j++)
        {
            double temp = 0.0;

            for (unsigned k = 0; k < m_colSize; k++)
            {
                temp += m_matrix[i][k] * B(k, j);
            }

            multip(i, j) = temp;
        }
    }

    return multip;
}

// Scalar Addition
Matrix Matrix::operator+(double scalar)
{
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++)
    {
        for (j = 0; j < m_colSize; j++)
        {
            result(i, j) = this->m_matrix[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
Matrix Matrix::operator-(double scalar)
{
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++)
    {
        for (j = 0; j < m_colSize; j++)
        {
            result(i, j) = this->m_matrix[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
Matrix Matrix::operator*(double scalar)
{
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++)
    {
        for (j = 0; j < m_colSize; j++)
        {
            result(i, j) = this->m_matrix[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
Matrix Matrix::operator/(double scalar)
{
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++)
    {
        for (j = 0; j < m_colSize; j++)
        {
            result(i, j) = this->m_matrix[i][j] / scalar;
        }
    }
    return result;
}

// Returns a mutable element when asked in the form A(x,y)
double &Matrix::operator()(const unsigned &rowNo, const unsigned &colNo)
{
    return m_matrix[rowNo][colNo];
}

// Returns a read-only element for const matrices.
const double &Matrix::operator()(const unsigned &rowNo, const unsigned &colNo) const
{
    return m_matrix[rowNo][colNo];
}

// No brainer - returns row #
unsigned Matrix::getRows() const
{
    return this->m_rowSize;
}

// returns col #
unsigned Matrix::getCols() const
{
    return this->m_colSize;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose()
{
    Matrix Transpose(m_colSize, m_rowSize, 0.0);
    for (unsigned i = 0; i < m_colSize; i++)
    {
        for (unsigned j = 0; j < m_rowSize; j++)
        {
            Transpose(i, j) = this->m_matrix[j][i];
        }
    }
    return Transpose;
}

//Setting proper width for each entry
void Matrix::prettyPrint(int precision, double eps) const
{
    std::cout << "Matrix:\n";

    int width = precision + 7;
    // 7 = sign + decimal point + up to 5 digits before decimal

    for (unsigned i = 0; i < m_rowSize; i++)
    {
        for (unsigned j = 0; j < m_colSize; j++)
        {
            double value = m_matrix[i][j];

            if (std::abs(value) < eps)
                value = 0.0;

            std::cout << "["
                      << std::setw(width)
                      << std::fixed
                      << std::setprecision(precision)
                      << value
                      << "]";
        }

        std::cout << '\n';
    }
}

//Approximating values of matrix
Matrix Matrix::clean(double eps) const
{
    Matrix result(*this);

    for (unsigned i = 0; i < m_rowSize; i++)
    {
        for (unsigned j = 0; j < m_colSize; j++)
        {
            if (std::abs(result(i, j)) < eps)
                result(i, j) = 0.0;
        }
    }

    return result;
}

// Prints the matrix beautifully
void Matrix::print() const
{
    cout << "Matrix: " << endl;
    for (unsigned i = 0; i < m_rowSize; i++)
    {
        for (unsigned j = 0; j < m_colSize; j++)
        {
            cout << "[" << m_matrix[i][j] << "] ";
        }
        cout << endl;
    }
}
