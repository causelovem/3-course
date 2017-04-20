#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <mpi.h>
#include <omp.h>
#include <complex>

using namespace std;
 
typedef std::complex <double> complexd;

int main(int argc, char *argv[])
{
    /*<fileA> <fileB> <e> <num of threads>*/
    if (argc != 5)
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

    //omp_set_num_threads(atoi(argv[4]));

    MPI_File fileA;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);

    MPI_File fileB;
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fileB);

    int n = 0, size = 0;
    double e = atof(argv[3]);

    MPI_File_seek(fileA, 0, MPI_SEEK_SET);
    MPI_File_read(fileA, &n, 1, MPI_INT, &status);

    size = (1 << n) / nProc;

    complexd *vecA = NULL, *vecB = NULL, *vecC = NULL;
    vecA = new complexd [size];
    vecB = new complexd [size];
    vecC = new complexd [size];

    int offset = 2 * (myRank * size) * sizeof(double) + sizeof(int);

    MPI_File_seek(fileA, offset, MPI_SEEK_SET);
    MPI_File_read(fileA, vecA, size, MPI_DOUBLE_COMPLEX, &status);

    double He[4];

    if (myRank == 0)
    {
        srand(time(0));

        double s = 0;
        for (int i = 0; i < 12; i++)
            s += (double)rand() / RAND_MAX;
        s = (s - 6) * e;

        He[0] = (1 / sqrt(2)) * (cos(s) - sin(s));
        He[1] = (1 / sqrt(2)) * (sin(s) + cos(s));
        He[2] = (1 / sqrt(2)) * (cos(s) + sin(s));
        He[3] = (1 / sqrt(2)) * (sin(s) - cos(s));
    }

    MPI_Bcast(He, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*
                      (1  1)
    H = (1 / sqrt(2)) (1 -1)

    */

    int ofs = myRank * size;

    float time = (float)clock();

    for (int k = 0; k < n; k++)
    {
        int mask = 1 << (n - (k + 1));
        int xor_ofs = ofs ^ mask;

        if (((xor_ofs < ofs + size) && (xor_ofs > ofs)) || (nProc == 1))
        {
            int size1 = size / 2;

            #pragma omp parallel for
            for (int i = 0; i < size1; i++)
            {
                int ormask = i ^ mask;

                complexd tmp1 = He[0] * vecA[i] + He[1] * vecA[ormask];
                complexd tmp2 = He[2] * vecA[i] + He[3] * vecA[ormask];
                complexd tmp1_clear = (1 / sqrt(2)) * (vecA[i] + vecA[ormask]);
                complexd tmp2_clear = (1 / sqrt(2)) * (vecA[i] - vecA[ormask]);

                vecA[i] = tmp1;
                vecA[ormask] = tmp2;
                vecC[i] = tmp1_clear;
                vecC[ormask] = tmp2_clear;
            }
        }
        else
        {
            int pairRank = xor_ofs / size;

            MPI_Sendrecv(vecA, size, MPI_DOUBLE_COMPLEX, pairRank, 0, vecB, size, MPI_DOUBLE_COMPLEX,
            pairRank, 0, MPI_COMM_WORLD, &status);

            if (myRank > pairRank)
            {
                #pragma omp parallel for
                for (int i = 0; i < size; i++)
                {
                    complexd tmp2 = He[2] * vecB[i] + He[3] * vecA[i];
                    complexd tmp2_clear = (1 / sqrt(2)) * (vecB[i] - vecA[i]);

                    vecA[i] = tmp2;
                    vecC[i] = tmp2_clear;
                }
            }
            else
            {
                #pragma omp parallel for
                for (int i = 0; i < size; i++)
                {
                    complexd tmp1 = He[0] * vecA[i] + He[1] * vecB[i];
                    complexd tmp1_clear = (1 / sqrt(2)) * (vecA[i] + vecB[i]);

                    vecA[i] = tmp1;
                    vecC[i] = tmp1_clear;
                }
            }
        }
    }

    float fin_time = (float)clock() - time;

    MPI_Reduce(&fin_time, &time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        cout << ">Time of computation = " << time / CLOCKS_PER_SEC << endl;
        cout << "N = " << n << endl;
        cout << ">Number of proc = " << nProc << endl;
        cout << ">Number of threads = " << atoi(argv[4]) << endl;
        cout << ">E = " << e << endl;
    }

    complexd tmp_sum = 0, sum = 0;

    /*for (int i = 0; i < size; i++)
        tmp_sum += norm(vecA[i]);

    MPI_Reduce(&tmp_sum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        sum = sqrt(sum);
        cout << sum << endl;
    }

    tmp_sum = sum = 0;*/

    //#pragma omp parallel for
    for (int i = 0; i < size; i++)
        tmp_sum += vecA[i] * conj(vecC[i]);

    MPI_Reduce(&tmp_sum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        complexd tmp (1, 0);
        sum = abs(sum);
        sum = tmp - sum * sum;
        cout << "1 - F = " << sum << endl;
    }

    
    /*if (myRank == 0)
    {
        MPI_File_seek(fileB, 0, MPI_SEEK_SET);
        MPI_File_write(fileB, &n, 1, MPI_INT, &status);
    }

    MPI_File_seek(fileB, offset, MPI_SEEK_SET);
    MPI_File_write(fileB, vecA, size, MPI_DOUBLE_COMPLEX, &status);*/

    MPI_File_close(&fileA);
    MPI_File_close(&fileB);

    delete [] vecA;
    delete [] vecB;
    delete [] vecC;

    MPI_Finalize();
    return 0;
}