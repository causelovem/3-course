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
    if (argc != 4)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    int left = atoi(argv[2]), right = atoi(argv[3]);

    if (right <= left)
    {
        cerr << ">Wrong arguments: left >= right." << endl;
        return -1;
    }

    MPI_Init(&argc, &argv);

    int nProc, myRank, root = 0, tag = 0;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    uint range[2];
    if (myRank == root)
    {
        ofstream outfile(argv[1]);
        if (outfile.is_open() == false)
        {
            cerr << ">Can not open outfile with such name." << endl;
            return -1;
        }

        uint res_num = 0;

        bool *primes = new bool [right + 1];

        primes[0] = primes[1] = 0;
        for (uint i = 2; i < right + 1; i++)
            primes[i] = 1;

        uint fin = sqrt(right) + 2;
        for (uint i = 2; i < fin; i++)
            if (primes[i] == 1)
                for (int j = i * i; j < right + 1; j += i)
                    primes[j] = 0;

        uint count = 0;
        for (uint i = 2; i < fin; i++)
            if (primes[i] == 1)
                count++;

        uint *tmp = new uint [count];
        uint cnt = 0;
        for (uint i = 2; i < fin; i++)
            if (primes[i] == 1)
                tmp[cnt++] = i;

        for (uint i = 0; i < count; i++)
        {
            outfile << tmp[i] << endl;
            res_num++;
        }

        MPI_Bcast(&count, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);

        MPI_Bcast(tmp, count, MPI_UNSIGNED, root, MPI_COMM_WORLD);

        int num = (right - fin ) / (nProc - 1);

        range[0] = fin - 1;
        range[1] = fin + num - 1;
        for (int i = 1; i < nProc; i++)
        {   
            MPI_Send(&range, 2, MPI_UNSIGNED, i, tag, MPI_COMM_WORLD);

            range[0] = range[1] + 1;
            range[1] += num + 1;
        }

        for (int i = 1; i < nProc; i++)
        {
            uint num = 0;
            MPI_Recv(&num, 1, MPI_UNSIGNED, i, tag, MPI_COMM_WORLD, &status);

            uint *primes_out = new uint [num];

            MPI_Recv(primes_out, num, MPI_UNSIGNED, i, tag, MPI_COMM_WORLD, &status);

            for (uint j = 0; j < num; j++)
            {
                outfile << primes_out[j] << endl;
                res_num++;
            }

            delete [] primes_out;
        }

        cout << ">Number of prime numbers = " << res_num << endl;

        delete [] primes;
        delete [] tmp;
        outfile.close();
    }
    else
    {
        uint count = 0;
        MPI_Bcast(&count, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
        uint *primes = new uint [count];

        MPI_Bcast(primes, count, MPI_UNSIGNED, root, MPI_COMM_WORLD);

        MPI_Recv(&range, 2, MPI_UNSIGNED, root, tag, MPI_COMM_WORLD, &status);
        uint cnt = range[1] - range[0] + 1;
        bool *tmp = new bool [cnt];

        for (uint i = 0; i < cnt; i++)
            tmp[i] = 1;

        for (uint i = 0; i < count; i++)
        {
            uint h = range[0] % primes[i];
            uint j = h == 0 ? 0 : primes[i] - h;
            for (; j < cnt; j += primes[i])
                tmp[j] = 0;
        }

        uint num = 0;
        for (uint i = 0; i < cnt; i++)
            if (tmp[i] == 1)
                num++;

        uint *primes_send = new uint [num];

        uint num_out = range[0];
        count = 0;
        for (uint i = 0; i < cnt; i++)
        {
            if (tmp[i] == 1)
                primes_send[count++] = num_out;
            num_out++;
        }

        MPI_Send(&num, 1, MPI_UNSIGNED, root, tag, MPI_COMM_WORLD);

        MPI_Send(primes_send, num, MPI_UNSIGNED, root, tag, MPI_COMM_WORLD);

        delete [] primes;
        delete [] tmp;
        delete [] primes_send;
    }

    MPI_Finalize();
    return 0;
}

    /*for (uint i = 2; i * i < right + 1; i++)
        if (primes[i] == 1)
            for (int j = i * i; j < right + 1; j += i)
                primes[j] = 0;*/