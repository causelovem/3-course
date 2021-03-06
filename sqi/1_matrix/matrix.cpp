#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

#include <vector>
#include <algorithm>

using namespace std;

template<typename T>
int ijk_multiple (std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    /*row * col*/    

    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    }

    std::vector<T> tmp_vec; 

    for (uint j = 0; j < matrixB[0].size(); j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < matrixA.size(); i++)
        matrixC.push_back(tmp_vec);

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    float time = 0;
    ofstream plot("plot.dat", ios_base::app);
    time = clock();
    for (uint i = 0; i < size1; i++)
    {
        for (uint j = 0; j < size2; j++)
        {
            T sum = 0.0;
            for (uint k = 0; k < size3; k++)
                sum += matrixA[i][k] * matrixB[k][j];
            matrixC[i][j] = sum;
        }
    }
    time = clock() - time;
    plot << 1 << " " << time / CLOCKS_PER_SEC << endl;

    plot.close();

    return 0;
}

template<typename T>
int ikj_multiple (std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    /*elem * row*/

    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    }

    std::vector<T> tmp_vec;

    for (uint j = 0; j < matrixB[0].size(); j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < matrixA.size(); i++)
        matrixC.push_back(tmp_vec);

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    float time = 0;
    ofstream plot("plot.dat", ios_base::app);
    time = clock();
    for (uint i = 0; i < size1; i++)
    {
        for (uint k = 0; k < size3; k++)
        {
            T elem = matrixA[i][k];
            for (uint j = 0; j < size2; j++)
                matrixC[i][j] += elem * matrixB[k][j];
        }
    }
    time = clock() - time;
    plot << 2 << " " << time / CLOCKS_PER_SEC << endl;

    plot.close();

    return 0;
}

template<typename T>
int jki_multiple (std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    /*elem * row*/

    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    }

    std::vector<T> tmp_vec;

    for (uint j = 0; j < matrixB[0].size(); j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < matrixA.size(); i++)
        matrixC.push_back(tmp_vec);

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    float time = 0;
    ofstream plot("plot.dat", ios_base::app);
    time = clock();
    for (uint j = 0; j < size2; j++)
    {
        for (uint k = 0; k < size3; k++)
        {
            T elem = matrixB[k][j];
            for (uint i = 0; i < size1; i++)
                matrixC[i][j] += elem * matrixA[i][k];
        }
    }
    time = clock() - time;
    plot << 3 << " " << time / CLOCKS_PER_SEC << endl;

    plot.close();

    return 0;
}

template<typename T>
int jik_multiple (std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    /*col * row*/

    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    }

    std::vector<T> tmp_vec;

    for (uint j = 0; j < matrixB[0].size(); j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < matrixA.size(); i++)
        matrixC.push_back(tmp_vec);

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    float time = 0;
    ofstream plot("plot.dat", ios_base::app);
    time = clock();
    for (uint j = 0; j < size2; j++)
    {
        for (uint i = 0; i < size1; i++)
        {
            T sum = 0.0;
            for (uint k = 0; k < size3; k++)
                sum += matrixA[i][k] * matrixB[k][j];
            matrixC[i][j] = sum;
        }
    }
    time = clock() - time;
    plot << 4 << " " << time / CLOCKS_PER_SEC << endl;

    plot.close();

    return 0;
}

template<typename T>
int kij_multiple (std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    /*elem * row*/

    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    }

    std::vector<T> tmp_vec;

    for (uint j = 0; j < matrixB[0].size(); j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < matrixA.size(); i++)
        matrixC.push_back(tmp_vec);

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    float time = 0;
    ofstream plot("plot.dat", ios_base::app);
    time = clock();
    for (uint k = 0; k < size3; k++)
    {
        for (uint i = 0; i < size1; i++)
        {
            T elem = matrixA[i][k];
            for (uint j = 0; j < size2; j++)
                matrixC[i][j] += elem * matrixB[k][j];
        }
    }
    time = clock() - time;
    plot << 5 << " " << time / CLOCKS_PER_SEC << endl;

    plot.close();

    return 0;
}

template<typename T>
int kji_multiple (std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    /*elem * row*/

    if (matrixA[0].size() != matrixB.size())
    {
        cerr << ">Can not multiple such matrix: wrong sizes." << endl;
        return -1;
    }

    std::vector<T> tmp_vec;

    for (uint j = 0; j < matrixB[0].size(); j++)
        tmp_vec.push_back(0);

    for (uint i = 0; i < matrixA.size(); i++)
        matrixC.push_back(tmp_vec);

    uint size1 = matrixA.size(), size2 = matrixB[0].size(), size3 = matrixB.size();

    float time = 0;
    ofstream plot("plot.dat", ios_base::app);
    time = clock();
    for (uint k = 0; k < size3; k++)
    {
        for (uint j = 0; j < size2; j++)
        {
            T elem = matrixB[k][j];
            for (uint i = 0; i < size1; i++)
                matrixC[i][j] += elem * matrixA[i][k];
        }
    }
    time = clock() - time;
    plot << 6 << " " << time / CLOCKS_PER_SEC << endl;

    plot.close();

    return 0;
}

template<typename T>
void read_matrix (std::vector < std::vector<T> > &matrix, ifstream &chanal)
{
    T num;
    uint row = 0, col = 0;
    std::vector<T> tmp_vec;

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

template<typename T>
int swich_type (int type, std::vector < std::vector<T> > &matrixA, std::vector < std::vector<T> > &matrixB, std::vector < std::vector<T> > &matrixC)
{
    if (type == 1)
    {
        if (ijk_multiple<T>(matrixA, matrixB, matrixC) == -1)
            return -1;
    }
    else
    if (type == 2)
    {
        if (ikj_multiple<T>(matrixA, matrixB, matrixC) == -1)
            return -1;
    }
    else
    if (type == 3)
    {
        if (jki_multiple<T>(matrixA, matrixB, matrixC) == -1)
            return -1;
    }
    else
    if (type == 4)
    {
        if (jik_multiple<T>(matrixA, matrixB, matrixC) == -1)
            return -1;
    }
    else
    if (type == 5)
    {
        if (kij_multiple<T>(matrixA, matrixB, matrixC) == -1)
            return -1;
    }
    else
    if (type == 6)
    {
        if (kji_multiple<T>(matrixA, matrixB, matrixC) == -1)
            return -1;
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    /*<A-file> <B-file> <C-file> <type>*/
    if (argc != 5)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    int ta = 0, tb = 0; /*0 - float, 1 - double*/
    int type = atoi(argv[4]); /*1 - ijk, 2 - ikj, 3 - jki, 4 - jik, 5 - kij, 6 - kji*/

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

    fileA >> ta;
    fileB >> tb;

    if ((ta == 0) && (tb == 0))
    {
        std::vector < std::vector<float> > matrixA;
        std::vector < std::vector<float> > matrixB;
        std::vector < std::vector<float> > matrixC;
        read_matrix(matrixA, fileA);
        read_matrix(matrixB, fileB);
        if (swich_type(type, matrixA, matrixB, matrixC) == -1)
        {
            fileA.close();
            fileB.close();
            fileC.close();
            return -1;
        }

        fileC << 0 << " " << matrixC.size() << " " << matrixC[0].size() << "\n";

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
    }
    else
    if (((ta == 1) && (tb == 1)) || ((ta == 1) && (tb == 0)) || ((ta == 0) && (tb == 1)))
    {
        std::vector < std::vector<double> > matrixA;
        std::vector < std::vector<double> > matrixB;
        std::vector < std::vector<double> > matrixC;
        read_matrix(matrixA, fileA);
        read_matrix(matrixB, fileB);
        if (swich_type(type, matrixA, matrixB, matrixC) == -1)
        {
            fileA.close();
            fileB.close();
            fileC.close();
            return -1;
        }

        fileC << 1 << " " << matrixC.size() << " " << matrixC[0].size() << "\n";

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
    }
    else
    {
        cerr << ">Unknown type of matrix." << endl;
        fileA.close();
        fileB.close();
        fileC.close();
        return -1;
    }

    fileA.close();
    fileB.close();
    fileC.close();
    return 0;
}