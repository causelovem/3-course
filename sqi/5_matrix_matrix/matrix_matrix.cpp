#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <mpi.h>

#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{
    /*<A-file> <B-file> <C-file>*/
    /*if (argc != 4)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }*/

    MPI_Init(&argc, &argv);

    MPI_File fileA;
    MPI_File fileB;
    MPI_File fileC;

    int nProc, myRank, root = 0, tag = 0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    double sq = pow(nProc,  1.0 / 3.0);
    int check_sq = trunc(sq);
    if ((check_sq = pow(check_sq, 3)) != nProc)
    {
        cout << ">Can not parallel program, because your proc number is not a 3-rd power of some N." << endl;

        MPI_Finalize();
        return -1;
    }
    int block_num = pow(nProc,  1.0 / 3.0);

    MPI_File_open(MPI_COMM_WORLD, "matrixB", MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);
    //MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);
    //MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileB);
    //MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY, MPI_INFO_NULL, &fileC);

    int size = 0;

    //MPI_File_set_view(file, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
    MPI_File_seek(fileA, 0, MPI_SEEK_SET);
    MPI_File_read(fileA, &size, 1, MPI_INT, &status);

    cout << myRank << " " << size << endl;

    MPI_File_close(&fileA);
    //MPI_File_close(&fileB);
    //MPI_File_close(&fileC);
    MPI_Finalize();
    return 0;
}