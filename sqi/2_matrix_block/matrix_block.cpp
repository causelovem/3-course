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

void read_matrix (float *&matrix, ifstream &chanal, uint &size)
{
    chanal >> size;
    chanal >> size;

    matrix = new float [size * size];

    for (uint i = 0; i < size * size; i++)
        chanal >> matrix[i];

    return;
}

static void do_block(uint n, uint M, uint N, uint K, float* A, float* B, float* C)
{
    for (uint i = 0; i < M; ++i)
    {
        const int iOffset = i*n;
        for (uint j = 0; j < N; ++j)
        {
            float cij = 0.0;
            for (uint k = 0; k < K; ++k)
                cij += A[iOffset+k] * B[k*n+j];
            C[iOffset+j] += cij;
        }
    }
}

void block_multiply(uint n , float* A, float* B, float* C, bool type, uint block_size)
{
    if (type == 1)
    {       
        for (uint i = 0; i < n; i+=block_size)
        {
            const uint iOffset = i * n;
            for (uint j = 0; j < n; j+=block_size)
            {
                for (uint k = 0; k < n; k+=block_size)
                {
                    uint M = min(block_size,n-i);
                    uint N = min(block_size,n-j);
                    uint K = min(block_size,n-k);

                    do_block(n, M, N, K, A + iOffset + k, B + k*n + j, C + iOffset + j);
                }
            }
        }
    }
    else
    {
        for (uint i = 0; i < n; i+=block_size)
        {
            const uint iOffset = i * n;
            for (uint k = 0; k < n; k+=block_size)
            {
                for (uint j = 0; j < n; j+=block_size)
                {
                    uint M = min(block_size,n-i);
                    uint N = min(block_size,n-j);
                    uint K = min(block_size,n-k);

                    do_block(n, M, N, K, A + iOffset + k, B + k*n + j, C + iOffset + j);
                }
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    /*<A-file> <B-file> <C-file> <block_size> <type>*/
    /*0 - ijk, 1 - ikj*/
    if (argc != 6)
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

    int block_size = atoi(argv[4]);
    int type = atoi(argv[5]);

    float *matrixA;
    float *matrixB;
    float *matrixC;

    uint size = 0;

    read_matrix(matrixA, fileA, size);

    read_matrix(matrixB, fileB, size);

    matrixC = new float [size * size];
    for (uint i = 0; i < size * size; i++)
        matrixC[i] = 0;

    int events[4] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_FP_OPS, PAPI_TOT_CYC};
    long_long values[4];

    if (PAPI_num_counters() < 4)
    {
        cout << ">No hardware counters here, or PAPI not supported." << endl;
        fileA.close();
        fileB.close();
        fileC.close();
        return -1;
    }

    if (PAPI_start_counters(events, 3) != PAPI_OK)
    {
        cout << ">PAPI faild to start counters." << endl;
        delete [] matrixA;
        delete [] matrixB;
        delete [] matrixC;
        fileA.close();
        fileB.close();
        fileC.close();
        return -1;
    }

    block_multiply(size, matrixA, matrixB, matrixC, type, block_size);

    if (PAPI_read_counters(values, 4) != PAPI_OK)
    {
        cout << ">PAPI faild to read counters." << endl;
        delete [] matrixA;
        delete [] matrixB;
        delete [] matrixC;
        fileA.close();
        fileB.close();
        fileC.close();
        return -1;
    }

    cout << ">Total hardware flops = " << (float)values[2] << endl;
    cout << ">L1 data cache misses is = " << (float)values[0] << endl;
    cout << ">L2 data cache misses is = " << (float)values[1] << endl;
    cout << ">Total cycles = " << (float)values[3] << endl;

    if (PAPI_stop_counters(values, 4) != PAPI_OK)
    {
        cout << ">PAPI faild to stop counters." << endl;
        delete [] matrixA;
        delete [] matrixB;
        delete [] matrixC;
        fileA.close();
        fileB.close();
        fileC.close();
        return -1;
    }

    fileC << size << " " << size << "\n";

    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            fileC << matrixC[i * size + j];
            if (j != size - 1)
                fileC << " ";
        }
        fileC << "\n";
    }

    delete [] matrixA;
    delete [] matrixB;
    delete [] matrixC;

    fileA.close();
    fileB.close();
    fileC.close();
    return 0;
}