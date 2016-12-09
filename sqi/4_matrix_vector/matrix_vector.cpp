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

    int nProc, myRank, root = 0, tag = 0;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == root)
    {
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

        cout << ">Number of proc = " << nProc << endl;

        std::vector < std::vector<double> > matrixA;

        double num;
        uint row = 0, col = 0;
        std::vector < double > tmp_vec;

        fileA >> row;
        fileA >> col;
        cout << ">MatrixA size = " << row << "x" << col << endl;
        for (uint i = 0; i < row; i++)
        {
            for (uint j = 0; j < col; j++)
            {
                fileA >> num;
                tmp_vec.push_back(num);
            }
            matrixA.push_back(tmp_vec);
            tmp_vec.clear();
        }

        std::vector < double > result(row);
        uint tmp_col = col;

        fileB >> col;
        std::vector < double > matrixB(col);
        cout << ">MatrixB size = " << col << endl;

        for (uint i = 0; i < col; i++)
            fileB >> matrixB[i];

        if (col != tmp_col)
        {
            cerr << ">Matrix and vector sizes are not comply with." << endl;
            fileA.close();
            fileB.close();
            fileC.close();
            MPI_Finalize();
            return -1;
        }

        if (nProc == 1)
        {
            float time = (float)clock();

            int res_count = 0;
            for (uint i = 0; i < row; i++)
            {
                double res = 0;
                for (uint j = 0; j < col; j++)
                    res += matrixA[i][j] * matrixB[j];
                result[res_count++] = res;
            }

            cout << ">Time of computation = " << ((float)clock() - time) / CLOCKS_PER_SEC << endl;

            fileC << row << "\n";
            for (uint i = 0; i < row; i++)
                fileC << result[i] << " ";

            fileA.close();
            fileB.close();
            fileC.close();
            MPI_Finalize();
            return 0;
        }

        if (nProc > row)
            nProc = row;

        float time = (float)clock();

        int div = row / nProc;
        int check = row % nProc;
        int count = div;

        double *vectorB = new double [col];

        for (uint i = 0; i < col; i++)
            vectorB[i] = matrixB[i];

        MPI_Bcast(&nProc, 1, MPI_INT, root, MPI_COMM_WORLD);

        MPI_Bcast(&col, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
        MPI_Bcast(vectorB, col, MPI_DOUBLE, root, MPI_COMM_WORLD);

        for (uint i = 1; i < nProc - 1; i++)
            MPI_Send(&div, 1, MPI_INT, i, tag, MPI_COMM_WORLD);

        int send_val = div + check;
        MPI_Send(&send_val, 1, MPI_INT, nProc - 1, tag, MPI_COMM_WORLD);

        for (uint i = 1; i < nProc - 1; i++)
        {
            for (uint j = count; j < count + div; j++)
            {
                for (uint k = 0; k < col; k++)
                    vectorB[k] = matrixA[j][k];
                MPI_Send(vectorB, col, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            }

            count += div;
        }

        for (uint i = count; i < count + send_val; i++)
        {
            for (uint k = 0; k < col; k++)
                vectorB[k] = matrixA[i][k];
            MPI_Send(vectorB, col, MPI_DOUBLE, nProc - 1, tag, MPI_COMM_WORLD);
        }

        int res_count = 0;
        for (uint i = 0; i < div; i++)
        {
            double res = 0;
            for (uint j = 0; j < col; j++)
                res += matrixA[i][j] * matrixB[j];
            result[res_count++] = res;
        }

        double *tmp = new double [div];

        for (uint i = 1; i < nProc - 1; i++)
        {
            MPI_Recv(tmp, div, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
            for (uint j = 0; j < div; j++)
                result[res_count++] = tmp[j];
        }

        double *send_vec = new double [send_val];
        MPI_Recv(send_vec, send_val, MPI_DOUBLE, nProc - 1, tag, MPI_COMM_WORLD, &status);

        for (uint j = 0; j < send_val; j++)
            result[res_count++] = send_vec[j];

        cout << ">Time of computation = " << ((float)clock() - time) / CLOCKS_PER_SEC << endl;

        fileC << row << "\n";
        for (uint i = 0; i < row; i++)
            fileC << result[i] << " ";

        fileA.close();
        fileB.close();
        fileC.close();

        delete [] vectorB;
        delete [] tmp;
        delete [] send_vec;
    }
    else
    {
        int proc_new;
        MPI_Bcast(&proc_new, 1, MPI_INT, root, MPI_COMM_WORLD);

        if (proc_new - 1 >= myRank)
        {
            uint col;
            MPI_Bcast(&col, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
            double *matrixB = new double [col];
            MPI_Bcast(matrixB, col, MPI_DOUBLE, root, MPI_COMM_WORLD);

            int div;
            MPI_Recv(&div, 1, MPI_INT, root, tag, MPI_COMM_WORLD, &status);

            double *res = new double [div];
            double **tmp = new double* [div];

            for (uint i = 0; i < div; i++)
            {
                tmp[i] = new double [col];
                MPI_Recv(tmp[i], col, MPI_DOUBLE, root, tag, MPI_COMM_WORLD, &status);
            }

            for (uint i = 0; i < div; i++)
                res[i] = 0;

            for (uint i = 0; i < div; i++)
                for (uint j = 0; j < col; j++)
                    res[i] += tmp[i][j] * matrixB[j];

            MPI_Send(res, div, MPI_DOUBLE, root, tag, MPI_COMM_WORLD);

            delete [] matrixB;
            delete [] res;

            for (uint i = 0; i < div; i++)
                delete [] tmp[i];

            delete [] tmp;
        }
    }
    
    MPI_Finalize();
    return 0;
}