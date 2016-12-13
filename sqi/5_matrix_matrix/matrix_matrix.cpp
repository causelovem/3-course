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
    if (argc != 4)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    MPI_Init(&argc, &argv);

    MPI_File fileA;
    MPI_File fileB;
    MPI_File fileC;

    int nProc, myRank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    double tmp = cbrtf(nProc);
    int dim = tmp;
    if (tmp - dim > 0.0000001)
    {
        if (myRank == 0)
            cout << ">Can not parallel program, because your proc number is not a 3-rd power of some N." << endl;

        MPI_Finalize();
        return -1;
    }

    int k = myRank / (dim * dim),
        j = (myRank % (dim * dim)) / dim,
        i = myRank - k * (dim * dim) - j * dim;


    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileA);
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fileB);
    MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fileC);

    float timeR = (float)clock();

    int size = 0;
    MPI_File_seek(fileA, 0, MPI_SEEK_SET);
    MPI_File_read(fileA, &size, 1, MPI_INT, &status);

    int nb_size = size / dim, mod = size % dim;
    
    MPI_Datatype matrix_block;
    MPI_Datatype matrix_block_R;
    MPI_Datatype matrix_block_B;
    MPI_Datatype matrix_block_RB;
    MPI_Type_vector(nb_size, nb_size, size, MPI_DOUBLE, &matrix_block);
    MPI_Type_vector(nb_size, nb_size + mod, size, MPI_DOUBLE, &matrix_block_R);
    MPI_Type_vector(nb_size + mod, nb_size, size, MPI_DOUBLE, &matrix_block_B);
    MPI_Type_vector(nb_size + mod, nb_size + mod, size, MPI_DOUBLE, &matrix_block_RB);

    MPI_Type_commit(&matrix_block);
    MPI_Type_commit(&matrix_block_R);
    MPI_Type_commit(&matrix_block_B);
    MPI_Type_commit(&matrix_block_RB);

    double *blockA = NULL, *blockB = NULL, *blockC = NULL;
    int sizeA = 0, sizeB = 0, sizeC = 0, rowA = 0, colA = 0, rowB = 0, colB = 0;

    int offsetA = ((i * nb_size) * size + (k * nb_size)) * sizeof(double) + sizeof(int);
    int offsetB = ((k * nb_size) * size + (j * nb_size)) * sizeof(double) + sizeof(int);
    int offsetC = ((i * nb_size) * size + (j * nb_size)) * sizeof(double) + sizeof(int);


    /*MATRIX A*/
    if ((i == dim - 1) && (k == dim - 1))
    {
        rowA = nb_size + mod;
        colA = nb_size + mod;
        sizeA = (nb_size + mod) * (nb_size + mod);
        blockA = new double [sizeA];
        MPI_File_set_view(fileA, offsetA, MPI_DOUBLE, matrix_block_RB, "native", MPI_INFO_NULL);
        sizeC = nb_size + mod;
    }
    else
    if ((i == dim - 1) && (k != dim - 1))
    {
        rowA = nb_size + mod;
        colA = nb_size;
        sizeA = (nb_size + mod) * nb_size;
        blockA = new double [sizeA];
        MPI_File_set_view(fileA, offsetA, MPI_DOUBLE, matrix_block_B, "native", MPI_INFO_NULL);
        sizeC = nb_size + mod;
    }
    else
    if ((i != dim - 1) && (k == dim - 1))
    {
        rowA = nb_size;
        colA = nb_size + mod;
        sizeA = nb_size * (nb_size + mod);
        blockA = new double [sizeA];
        MPI_File_set_view(fileA, offsetA, MPI_DOUBLE, matrix_block_R, "native", MPI_INFO_NULL);
        sizeC = nb_size;
    }
    else
    {
        rowA = nb_size;
        colA = nb_size;
        sizeA = nb_size * nb_size;
        blockA = new double [sizeA];
        MPI_File_set_view(fileA, offsetA, MPI_DOUBLE, matrix_block, "native", MPI_INFO_NULL);
        sizeC = nb_size;
    }
    MPI_File_read(fileA, blockA, sizeA, MPI_DOUBLE, &status);


    /*MATRIX B*/
    if ((k == dim - 1) && (j == dim - 1))
    {
        rowB = nb_size + mod;
        colB = nb_size + mod;
        sizeB = (nb_size + mod) * (nb_size + mod);
        blockB = new double [sizeB];
        MPI_File_set_view(fileB, offsetB, MPI_DOUBLE, matrix_block_RB, "native", MPI_INFO_NULL);
        sizeC *= nb_size + mod;
    }
    else
    if ((k == dim - 1) && (j != dim - 1))
    {
        rowB = nb_size + mod;
        colB = nb_size;
        sizeB = (nb_size + mod) * nb_size;
        blockB = new double [sizeB];
        MPI_File_set_view(fileB, offsetB, MPI_DOUBLE, matrix_block_B, "native", MPI_INFO_NULL);
        sizeC *= nb_size;
    }
    else
    if ((k != dim - 1) && (j == dim - 1))
    {
        rowB = nb_size;
        colB = nb_size + mod;
        sizeB = nb_size * (nb_size + mod);
        blockB = new double [sizeB];
        MPI_File_set_view(fileB, offsetB, MPI_DOUBLE, matrix_block_R, "native", MPI_INFO_NULL);
        sizeC *= nb_size + mod;
    }
    else
    {
        rowB = nb_size;
        colB = nb_size;
        sizeB = nb_size * nb_size;
        blockB = new double [sizeB];
        MPI_File_set_view(fileB, offsetB, MPI_DOUBLE, matrix_block, "native", MPI_INFO_NULL);
        sizeC *= nb_size;
    }
    MPI_File_read(fileB, blockB, sizeB, MPI_DOUBLE, &status);


    blockC = new double [sizeC];
    double *blockR = new double [sizeC];

    for (int l = 0; l < sizeC; l++)
        blockC[l] = 0;

    float timeC = (float)clock();

    for (int l = 0; l < rowA; l++)
        for (int t = 0; t < colB; t++)
            for (int s = 0; s < colA; s++)
                blockC[l * colB + t] += blockA[l * colA + s] * blockB[t + s * colB];

    timeC = ((float)clock() - timeC) / CLOCKS_PER_SEC;

    int redRank = myRank % (dim * dim);
    int process_ranks[dim];

    for (int l = 0; l < dim; l++)
        process_ranks[l] = redRank + dim * dim * l;

    int rank_in_group;
    MPI_Group group_world, sub_group;
    MPI_Comm sub_group_comm;
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, dim, process_ranks, &sub_group);
    MPI_Comm_create(MPI_COMM_WORLD, sub_group, &sub_group_comm);
    MPI_Comm_rank(sub_group_comm, &rank_in_group);
    MPI_Reduce(blockC, blockR, sizeC, MPI_DOUBLE, MPI_SUM, 0, sub_group_comm);


    /*MATRIX C*/
    if (myRank == 0)
    {
        MPI_File_seek(fileC, 0, MPI_SEEK_SET);
        MPI_File_write(fileC, &size, 1, MPI_INT, &status);
    }

    if ((i == dim - 1) && (j == dim - 1))
    {
        MPI_File_set_view(fileC, offsetC, MPI_DOUBLE, matrix_block_RB, "native", MPI_INFO_NULL);
    }
    else
    if ((i == dim - 1) && (j != dim - 1))
    {
        MPI_File_set_view(fileC, offsetC, MPI_DOUBLE, matrix_block_B, "native", MPI_INFO_NULL);
    }
    else
    if ((i != dim - 1) && (j == dim - 1))
    {
        MPI_File_set_view(fileC, offsetC, MPI_DOUBLE, matrix_block_R, "native", MPI_INFO_NULL);
    }
    else
    {
        MPI_File_set_view(fileC, offsetC, MPI_DOUBLE, matrix_block, "native", MPI_INFO_NULL);
    }    

    if (myRank < dim * dim)
        MPI_File_write(fileC, blockR, sizeC, MPI_DOUBLE, &status);

    timeR = ((float)clock() - timeR) / CLOCKS_PER_SEC;

    float redtimeC = 0, redtimeR = 0;

    MPI_Reduce(&timeC, &redtimeC, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeR, &redtimeR, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        cout << ">Matrix size = " << size << endl;
        cout << ">DIM = " << dim << endl;
        cout << ">Time of computation = " << redtimeC << endl;
        cout << ">Programm time = " << redtimeR << endl;
    }

    delete [] blockA;
    delete [] blockB;
    delete [] blockC;
    delete [] blockR;
    blockA = blockB = blockC = blockR = NULL;

    MPI_File_close(&fileA);
    MPI_File_close(&fileB);
    MPI_File_close(&fileC);
    MPI_Finalize();
    return 0;
}