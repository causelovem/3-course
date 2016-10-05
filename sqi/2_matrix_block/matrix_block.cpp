#include <papi.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

#include <vector>
#include <algorithm>

using namespace std;

void read_matrix (std::vector < std::vector<float> > &matrix, ifstream &chanal)
{
    float num;
    uint row = 0, col = 0;
    std::vector<float> tmp_vec;

    chanal >> row;
    chanal >> col;
    for (uint i = 0; i < row; i++)
    {
        for (uint j = 0; j < col; j++)
        {
            chanal >> num;
            tmp_vec.push_back(num);
        }
        matrix.push_back(tmp_vec);
        tmp_vec.clear();
    }

    return;
}

int block_multiple (std::vector < std::vector<float> > &matrixA, std::vector < std::vector<float> > &matrixB, std::vector < std::vector<float> > &matrixC, uint block_size)
{
    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    } 

    std::vector<float> tmp_vec;

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    for (uint j = 0; j < size2; j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < size1; i++)
        matrixC.push_back(tmp_vec);

    for (uint i = 0; i < size1; i = i + block_size)
        for (uint j = 0; j < size2; j = j + block_size)
            for (uint k = 0; k < size3; k = k + block_size)
            {
                uint l = i + block_size, t = j + block_size, n = k + block_size;
                for (uint i1 = 0; i1 < l; i1++)
                {
                    for (uint j1 = 0; j1 < t; j1++)
                    {
                        float elem = matrixA[i1][j1];
                        for (uint k1 = 0; k1 < n; k1++)
                            matrixC[i1][k1] += elem * matrixB[j1][k1];
                    }
                }
            }

    return 0;
}

int main(int argc, char const *argv[])
{
    /*<A-file> <B-file> <C-file> <block size>*/
    if (argc != 5)
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

    ofstream fileC(argv[3]);
    if (fileC.is_open() == false)
    {
        cerr << ">Can not open C-file with such name." << endl;
        fileA.close();
        fileB.close();
        return -1;
    }

    uint block_size = atoi(argv[4]);

    std::vector < std::vector<float> > matrixA;
    std::vector < std::vector<float> > matrixB;
    std::vector < std::vector<float> > matrixC;
    read_matrix(matrixA, fileA);
    read_matrix(matrixB, fileB);

    if (block_multiple(matrixA, matrixB, matrixC, block_size) == -1)
    {
        fileA.close();
        fileB.close();
        fileC.close();
        return -1;
    }

    fileC << matrixC.size() << " " << matrixC[0].size() << "\n";

    for (uint i = 0; i < matrixC.size(); i++)
    {
        for (uint j = 0; j < matrixC[0].size(); j++)
        {
            fileC << matrixC[i][j];
            if (j != matrixC[0].size() - 1)
                fileC << " ";
        }
        fileC << "\n";
    }

    fileA.close();
    fileB.close();
    fileC.close();
    return 0;
}