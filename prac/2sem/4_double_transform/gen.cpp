#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include <mpi.h>
#include <omp.h>
#include <complex>

using namespace std;

typedef std::complex <double> complexd;


int main(int argc, char *argv[])
{
    /*<A-file> <N>*/
    if (argc != 3)
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

    int n = atoi(argv[2]), size = 1 << n;

    int size_tmp = size / nProc;

    MPI_File fileA;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fileA);

    int rand_val = 0;

    if (myRank == 0)
    {
        srand(time(0));

        for (uint i = 1; i < nProc; i++)
        {
            rand_val = rand();
            MPI_Send(&rand_val, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&rand_val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        srand(rand_val);
    }

    complexd *tmp_vec = new complexd [size_tmp];

    complexd tmp_sum = 0, sum = 0;
    for (int i = 0; i < size_tmp; i++)
    {
        double re = (double)(rand()) / RAND_MAX * 20 - 10;
        double im = (double)(rand()) / RAND_MAX * 20 - 10;

        complexd c(re, im);

        tmp_vec[i] = c;

        //tmp_sum += c * c;
        tmp_sum += norm(c);
    }

    MPI_Reduce(&tmp_sum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0)
        sum = sqrt(sum);

    MPI_Bcast(&sum, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    #pragma omp parallel for
    for (int i = 0; i < size_tmp; i++)
        tmp_vec[i] /= sum;

    /*tmp_sum = sum = 0;
    for (int i = 0; i < size_tmp; i++)
        tmp_sum += norm(tmp_vec[i]);

    MPI_Reduce(&tmp_sum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        sum = sqrt(sum);
        cout << sum << endl;
    }*/

    if (myRank == 0)
    {
        MPI_File_seek(fileA, 0, MPI_SEEK_SET);
        MPI_File_write(fileA, &n, 1, MPI_INT, &status);
    }

    MPI_File_seek(fileA, 2 * (myRank * size_tmp) * sizeof(double) + sizeof(int), MPI_SEEK_SET);
    MPI_File_write(fileA, tmp_vec, size_tmp, MPI_DOUBLE_COMPLEX, &status);

    delete [] tmp_vec;

    MPI_File_close(&fileA);

    MPI_Finalize();
    return 0;
}