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
    /*<fileA> <fileB> <k1> <k2>*/
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

    MPI_File fileA;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);

    MPI_File fileB;
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fileB);

    int n = 0, k1 = atoi(argv[3]), k2 = atoi(argv[4]), size = 0;

    MPI_File_seek(fileA, 0, MPI_SEEK_SET);
    MPI_File_read(fileA, &n, 1, MPI_INT, &status);

    if ((k1 > n) || (k2 > n))
    {
        cerr << ">Number of bit is too big (K > N)." << endl;
        MPI_File_close(&fileA);
        MPI_File_close(&fileB);
        return -1;
    }

    size = (1 << n) / nProc;
    int full_size = (1 << n);

    complexd *vecA = NULL, *res = NULL;
    vecA = new complexd [full_size];
    res = new complexd [size];

    int offset = 2 * (myRank * size) * sizeof(double) + sizeof(int), mask1 = 1 << (n - k1), mask2 = 1 << (n - k2);
    //int ofs = myRank * size, xor_ofs1 = ofs ^ mask1, xor_ofs2 = ofs ^ mask2, xor_ofs3 = ofs ^ (mask1 | mask2);
    int ofs = myRank * size;


    MPI_File_seek(fileA, sizeof(int), MPI_SEEK_SET);

    MPI_File_read(fileA, vecA, full_size, MPI_DOUBLE_COMPLEX, &status);

    //MPI_File_seek(fileA, offset, MPI_SEEK_SET);

    //MPI_File_read(fileA, vecA, size, MPI_DOUBLE_COMPLEX, &status);

    complexd matrix[4][4] = {{1, 0, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 0, 1},
                             {0, 0, 1, 0}};

    //double time = MPI_Wtime();
    double time = clock();                         
    //cout << time << " " << clock() << endl;
    #pragma omp parallel for
    for (int i = ofs; i < ofs + size; i++)
    {
        int i00 = (i & ~mask1) & ~mask2;
        int i01 = (i & ~mask1) | mask2;
        int i10 = (i | mask1) & ~mask2;
        int i11 = (i | mask1) | mask2;

        int line = (((i & mask1) >> (n - k1)) << 1) + ((i & mask2) >> (n - k2));

        res[i - ofs] = matrix[line][0] * vecA[i00] + matrix[line][1] * vecA[i01]
        + matrix[line][2] * vecA[i10] + matrix[line][3] * vecA[i11];
    }


    //cout << (float)MPI_Wtime() << " " << clock() << endl;
    //double fin_time = MPI_Wtime() - time;
    double fin_time = clock() - time;

    MPI_Reduce(&fin_time, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        //cout << ">Time of computation = " << time << endl;
        cout << ">Time of computation = " << time / CLOCKS_PER_SEC << endl;
        cout << "N = " << (1 << n) << "\nK1 = " << k1 << "\nK2 = " << k2 << endl;
        cout << ">Number of proc = " << nProc << endl;
    }

    /*if (myRank == 0)
    {
        MPI_File_seek(fileB, 0, MPI_SEEK_SET);
        MPI_File_write(fileB, &n, 1, MPI_INT, &status);
    }

    MPI_File_seek(fileB, offset, MPI_SEEK_SET);
    MPI_File_write(fileB, res, size, MPI_DOUBLE_COMPLEX, &status);*/

    MPI_File_close(&fileA);
    MPI_File_close(&fileB);

    delete [] vecA;
    delete [] res;

    MPI_Finalize();
    return 0;
}



//float time = 0;
    /*if ((((xor_ofs1 < ofs + size) && (xor_ofs1 > ofs)) && ((xor_ofs2 < ofs + size) && (xor_ofs2 > ofs))) || (nProc == 1))
    {
        //cout << "!!!!" << endl;
        time = (float)clock();

        #pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            int i00 = (i & ~mask1) & ~mask2;
            int i01 = (i & ~mask1) | mask2;
            int i10 = (i | mask1) & ~mask2;
            int i11 = (i | mask1) | mask2;

            int line = (((i & mask1) >> (n - k1)) << 1) + ((i & mask2) >> (n - k2));

            res[i] = matrix[line][0] * vecA[i00] + matrix[line][1] * vecA[i01]
            + matrix[line][2] * vecA[i10] + matrix[line][3] * vecA[i11];
        }

        float fin_time = (float)clock() - time;

        MPI_Reduce(&fin_time, &time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    }
    else
    if (((xor_ofs1 < ofs + size) && (xor_ofs1 > ofs)) || ((xor_ofs2 < ofs + size) && (xor_ofs2 > ofs)))
    {
        //cout << "????" << endl;
        time = (float)clock();

        int pairRank1 = xor_ofs1 / size;
        int pairRank2 = xor_ofs2 / size;
        int pairRank = 0;

        complexd *tmp1 = NULL;
        tmp1 = new complexd [size];

        if (pairRank2 == myRank)
        {
            MPI_Sendrecv(vecA, size, MPI_DOUBLE_COMPLEX, pairRank1, 0, tmp1, size, MPI_DOUBLE_COMPLEX,
            pairRank1, 0, MPI_COMM_WORLD, &status);
            pairRank = pairRank1;
        }
        else
        {
            MPI_Sendrecv(vecA, size, MPI_DOUBLE_COMPLEX, pairRank2, 0, tmp1, size, MPI_DOUBLE_COMPLEX,
            pairRank2, 0, MPI_COMM_WORLD, &status);
            pairRank = pairRank2;
        }

        if (myRank < pairRank)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
            {
                int i00 = (i & ~mask1) & ~mask2;
                int i01 = (i & ~mask1) | mask2;
                //int i10 = (i | mask1) & ~mask2;
                //int i11 = (i | mask1) | mask2;
                
                int line = (((i & mask1) >> (n - k1)) << 1) + ((i & mask2) >> (n - k2));

                res[i] = matrix[line][0] * vecA[i00] + matrix[line][1] * vecA[i01]
                + matrix[line][2] * tmp1[i00] + matrix[line][3] * tmp1[i01];
            }
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
            {
                int i00 = (i & ~mask1) & ~mask2;
                int i01 = (i & ~mask1) | mask2;
                //int i10 = (i | mask1) & ~mask2;
                //int i11 = (i | mask1) | mask2;
                
                int line = (((i & mask1) >> (n - k1)) << 1) + ((i & mask2) >> (n - k2));

                res[i] = matrix[line][0] * tmp1[i00] + matrix[line][1] * tmp1[i01]
                + matrix[line][2] * vecA[i00] + matrix[line][3] * vecA[i01];
            }
        }

        float fin_time = (float)clock() - time;

        MPI_Reduce(&fin_time, &time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

        delete [] tmp1;
    }
    else
    {
        //cout << "&&&&" << myRank << endl;
        time = (float)clock();

        int pairRank1 = xor_ofs1 / size;
        int pairRank2 = xor_ofs2 / size;
        int pairRank3 = xor_ofs3 / size;

        complexd *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp = NULL;
        tmp1 = new complexd [size];
        tmp2 = new complexd [size];
        tmp3 = new complexd [size];

        int rank[3] = {pairRank1, pairRank2, pairRank3};

        for (int i = 0; i < 2; i++)
            for (int j = i + 1; j < 3; j++)
                if (rank[j] < rank[i])
                {
                    int tmp = rank[j];
                    rank[j] = rank[i];
                    rank[i] = tmp;
                }

        MPI_Sendrecv(vecA, size, MPI_DOUBLE_COMPLEX, rank[0], 0, tmp1, size, MPI_DOUBLE_COMPLEX,
        rank[0], 0, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(vecA, size, MPI_DOUBLE_COMPLEX, rank[1], 0, tmp2, size, MPI_DOUBLE_COMPLEX,
        rank[1], 0, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(vecA, size, MPI_DOUBLE_COMPLEX, rank[2], 0, tmp3, size, MPI_DOUBLE_COMPLEX,
        rank[2], 0, MPI_COMM_WORLD, &status);

        if ((myRank < rank[1]) && (myRank > rank[0]))
        {
            tmp = vecA;
            vecA = tmp1;
            tmp1 = tmp;
            tmp = NULL;
        }
        else
        if ((myRank < rank[2]) && (myRank > rank[1]))
        {
            tmp = vecA;
            vecA = tmp1;
            tmp1 = tmp2;
            tmp2 = tmp;
            tmp = NULL;
        }
        else
        if (myRank > rank[2])
        {
            tmp = vecA;
            vecA = tmp1;
            tmp1 = tmp2;
            tmp2 = tmp3;
            tmp3 = tmp;
            tmp = NULL;
        }

        #pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            int line = (((i & mask1) >> (n - k1)) << 1) + ((i & mask2) >> (n - k2));

            res[i] = matrix[line][0] * vecA[i] + matrix[line][1] * tmp1[i]
            + matrix[line][2] * tmp2[i] + matrix[line][3] * tmp3[i];

            cout << res[i] << endl;
        }

        float fin_time = (float)clock() - time;

        MPI_Reduce(&fin_time, &time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

        delete [] tmp1;
        delete [] tmp2;
        delete [] tmp3;
    }*/