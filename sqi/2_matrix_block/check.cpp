#include <octave/oct.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[])
{
    if (argc != 4)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    ifstream fileA(argv[1]);
    if (fileA.is_open() == false)
    {
        cerr << ">Can not open A-file with such name." << endl;
        return -1;
    }

    ifstream fileB(argv[2]);
    if (fileB.is_open() == false)
    {
        cerr << ">Can not open B-file with such name." << endl;
        fileA.close();
        return -1;
    }

    ifstream fileC(argv[3]);
    if (fileC.is_open() == false)
    {
        cerr << ">Can not open C-file with such name." << endl;
        fileA.close();
        fileB.close();
        return -1;
    }

    int row = 0, col = 0;

    fileA >> row;
    fileA >> col;

    Matrix A = Matrix(row, col);
    for(uint i = 0; i < row; i++)
        for (uint j = 0; j < col; j++)
            fileA >> A(i, j);

    fileB >> row;
    fileB >> col;

    Matrix B = Matrix(row, col);
    for(uint i = 0; i < row; i++)
        for (uint j = 0; j < col; j++)
            fileB >> B(i, j);

    fileC >> row;
    fileC >> col;

    Matrix C = Matrix(row, col);
    for(uint i = 0; i < row; i++)
        for (uint j = 0; j < col; j++)
            fileC >> C(i, j);

    Matrix res = A * B;
    double eps = 1;

    for(uint i = 0; i < row; i++)
        for (uint j = 0; j < col; j++)
            if(fabs(fabs(res(i, j)) - fabs(C(i, j))) > eps)
            {
                cerr << ">Matrixs are not equal." << endl;
                fileA.close();
                fileB.close();
                fileC.close();
                return -1;
            }
        /*{
            cout << fabs(res(i, j)) << " "  << fabs(C(i, j)) << endl;
        }*/

    cerr << ">Matrixs are equal*_*" << endl;

    fileA.close();
    fileB.close();
    fileC.close();
    return 0;
}