#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <mpi.h>
#include <complex>

using namespace std;

#define EPS 0.001
 
typedef std::complex <double> complexd;

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (((nProc != 0) && ((nProc & (~nProc + 1)) != nProc)))
    {
        if (myRank == 0)
            cout << ">Can not parallel program, because your proc number is not some 2 ^ N." << endl;

        MPI_Finalize();
        return -1;
    }

    MPI_File fileA;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);

    MPI_File fileB;
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileB);

    int n = 0, k = atoi(argv[3]), size = 0;

    MPI_File_seek(fileA, 0, MPI_SEEK_SET);
    MPI_File_read(fileA, &n, 1, MPI_INT, &status);

    if (k > (1 << n))
    {
        cerr << ">Number of bit is too big (K > N)." << endl;
        MPI_File_close(&fileA);
        MPI_File_close(&fileB);
        return -1;
    }

    size = (1 << n) / nProc;

    complexd *vecA = NULL, *vecB = NULL;
    vecA = new complexd [size];
    vecB = new complexd [size];

    int offset = 2 * (myRank * size) * sizeof(double) + sizeof(int), mask = 1 << (n - k);
    int ofs = myRank * size, xor_ofs = ofs ^ mask;

    MPI_File_seek(fileA, offset, MPI_SEEK_SET);
    MPI_File_seek(fileB, offset, MPI_SEEK_SET);

    MPI_File_read(fileA, vecA, size, MPI_DOUBLE_COMPLEX, &status);
    MPI_File_read(fileB, vecB, size, MPI_DOUBLE_COMPLEX, &status);

    /*
                      (1  1)
    H = (1 / sqrt(2)) (1 -1)

    */

    if (((xor_ofs < ofs + size) && (xor_ofs > ofs)) || (nProc == 1))
    {
        float time = (float)clock();

        int size1 = size / 2;
        for (int i = 0; i < size1; i++)
        {
            int ormask = i ^ mask;

            complexd tmp1 = (1 / sqrt(2)) * (vecB[i] + vecB[ormask]);
            complexd tmp2 = (1 / sqrt(2)) * (vecB[i] - vecB[ormask]);

            vecB[i] = tmp1;
            vecB[ormask] = tmp2;
        }

        float fin_time = (float)clock() - time;

        MPI_Reduce(&fin_time, &time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

        if (myRank == 0)
        {
            cout << ">Time of computation = " << time / CLOCKS_PER_SEC << endl;
            cout << "N = " << (1 << n) << "\nK = " << k << endl;
            cout << ">Number of proc = " << nProc << endl;
        }
    }
    else
    {
        float time = (float)clock();

        int pairRank = xor_ofs / size;

        complexd *vecC = NULL;
        vecC = new complexd [size];

        MPI_Sendrecv(vecB, size, MPI_DOUBLE_COMPLEX, pairRank, 0, vecC, size, MPI_DOUBLE_COMPLEX,
        pairRank, 0, MPI_COMM_WORLD, &status);

        if (myRank > pairRank)
        {
            for (int i = 0; i < size; i++)
            {
                complexd tmp2 = (1 / sqrt(2)) * (vecC[i] - vecB[i]);

                vecB[i] = tmp2;
            }
        }
        else
        {
            for (int i = 0; i < size; i++)
            {
                complexd tmp1 = (1 / sqrt(2)) * (vecB[i] + vecC[i]);

                vecB[i] = tmp1;
            }
        }

        float fin_time = (float)clock() - time;

        MPI_Reduce(&fin_time, &time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

        if (myRank == 0)
        {
            cout << ">Time of computation = " << time / CLOCKS_PER_SEC << endl;
            cout << "N = " << (1 << n) << "\nK = " << k << endl;
            cout << ">Number of proc = " << nProc << endl;
        }

        delete [] vecC;
    }

    int check = 0, fin_check = 0;
    for (int i = 0; i < size; i++)
        if ((abs(vecA[i]) - abs(vecB[i])) > EPS)
            check++;

    MPI_Reduce(&check, &fin_check, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        if (fin_check == 0)
            cout << "YES *_*\nvecA == vecB" << endl;
        else
            cout << "NO :(" << endl;
    }

    MPI_File_close(&fileA);
    MPI_File_close(&fileB);

    delete [] vecA;
    delete [] vecB;

    MPI_Finalize();
    return 0;
}