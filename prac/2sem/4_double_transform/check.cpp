#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <mpi.h>
#include <complex>

using namespace std;

#define EPS 0.1
 
typedef std::complex <double> complexd;

int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    MPI_File fileA;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);

    MPI_File fileB;
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileB);

    int n = 0, k1 = atoi(argv[4]), k2 = atoi(argv[5]), size = 0;

    MPI_File_seek(fileA, 0, MPI_SEEK_SET);
    MPI_File_read(fileA, &n, 1, MPI_INT, &status);

    MPI_File_seek(fileB, 0, MPI_SEEK_SET);
    MPI_File_read(fileB, &n, 1, MPI_INT, &status);

    size = (1 << n);

    complexd *vecA = NULL, *res = NULL, *vecB = NULL;
    vecA = new complexd [size];
    res = new complexd [size];
    vecB = new complexd [size];


    complexd matrix[4][4] = {{1, 0, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 0, 1},
                             {0, 0, 1, 0}};

    MPI_File_seek(fileA, sizeof(int), MPI_SEEK_SET);
    MPI_File_seek(fileB, sizeof(int), MPI_SEEK_SET);

    MPI_File_read(fileA, vecA, size, MPI_DOUBLE_COMPLEX, &status);
    MPI_File_read(fileB, res, size, MPI_DOUBLE_COMPLEX, &status);

    int mask1 = 1 << (n - k1), mask2 = 1 << (n - k2);

    for (int i = 0; i < size; i++)
    {
        int i00 = (i & ~mask1) & ~mask2;
        int i01 = (i & ~mask1) | mask2;
        int i10 = (i | mask1) & ~mask2;
        int i11 = (i | mask1) | mask2;

        int line = (((i & mask1) >> (n - k1)) << 1) + ((i & mask2) >> (n - k2));

        vecB[i] = matrix[line][0] * vecA[i00] + matrix[line][1] * vecA[i01]
        + matrix[line][2] * vecA[i10] + matrix[line][3] * vecA[i11];
    }

    int check = 0;
    for (int i = 0; i < size; i++)
    {
        if (abs(res[i]) - abs(vecB[i]) > EPS)
        {
            check = 1;
            break;
        }
    }

    if (check == 0)
        cout << "YES *_*\nvecA == vecB" << endl;
    else
        cout << "NO :(" << endl;

    MPI_File_close(&fileA);
    MPI_File_close(&fileB);

    delete [] vecA;
    delete [] vecB;
    delete [] res;

    MPI_Finalize();
    return 0;
}