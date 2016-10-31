#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <mpi.h>

#include <vector>
#include <algorithm>

#define GEN_BEGIN 300
#define GEN_END 1000
#define PERS 0.3
#define MUT_PERS 3
#define EXCHANGE 10

using namespace std;

struct item
{
    int vol;
    int val;
};

struct parents
{
    uint par1;
    uint par2;
};

struct generation
{
    uint gen;
    int gen_fit;
};

struct scomp
{
    bool operator () (generation fit1, generation fit2)
    {
        return fit1.gen_fit < fit2.gen_fit;
    }
} comp;

void mutation (uint &gen, const uint quan)
{
    uint pos = 0, point = 1, tmp_gen = 0;

    pos = (uint)(rand()) % quan + 1;

    tmp_gen = gen;

    point = point << pos;

    tmp_gen = tmp_gen & point;

    tmp_gen = tmp_gen ^ point;

    point = ~point;

    gen = gen & point;
    gen = gen | tmp_gen;

    return;
}

int fitness (const generation &individual, const vector< item > &items, const int volume)
{
    int fit = 0, v = 0;
    uint one = 1, tmp = 0, count = 0;
    uint quan = items.size();

    for (int i = quan - 1; i > -1; i--)
    {
        tmp = individual.gen & one;
        tmp = tmp >> count;
        fit +=  tmp * items[i].val;
        v += tmp * items[i].vol;
        one = one << 1;
        count++;
    }

    if (v > volume)
        fit = 0;

    return fit;
}

uint cross (parents par, const uint quan)
{
    uint point = 0, res = 0;
    uint top = 1;
    top = (top << quan) - 1;

    point = (uint)(rand()) % top + 1;

    par.par1 = par.par1 >> point;
    par.par2 = par.par2 << (sizeof(uint) * 8 - point);

    par.par1 = par.par1 << point;
    par.par2 = par.par2 >> (sizeof(uint) * 8 - point);

    res = par.par1 | par.par2;

    point = (uint)(rand()) % 100 + 1;
    if (point < MUT_PERS)
        mutation(res, quan);

    return res;
}

void selection (const vector< item > &items, vector<generation> &generat, const int volume, const int to, const int from)
{
    vector<parents> parent;
    uint quan = items.size();

    uint pair = GEN_BEGIN * PERS;

    int tag1 = 0, tag2 = 0;
    MPI_Status status;

    for (uint i = 0; i < GEN_END; i++)
    {
        for (uint j = 0; j < GEN_BEGIN; j++)
            generat[j].gen_fit = fitness(generat[j], items, volume);
        sort(generat.begin(), generat.end(), comp);

        parent.clear();
        parent.resize(pair);
        for (uint j = 0; j < pair; j++)
        {
            parent[j].par1 = 0;
            parent[j].par2 = 0;
        }

        for (uint j = 0; j < pair; j++)
        {
            parent[j].par1 = generat[(uint)(rand()) % GEN_BEGIN].gen;
            parent[j].par2 = generat[(uint)(rand()) % GEN_BEGIN].gen;

            while (parent[j].par1 == parent[j].par2)
                parent[j].par2 = generat[(uint)(rand()) % GEN_BEGIN].gen;
        }

        for (uint j = 0; j < pair; j++)
            generat[j].gen = cross(parent[j], quan);

        if ((i != 0) && ((i % EXCHANGE) == 0))
        {
            sort(generat.begin(), generat.end(), comp);

            uint *send = new uint [pair];
            uint *recv = new uint [pair];

            int size = generat.size();
            for (uint j = 0; j < pair; j++)
                send[j] = generat[size - 1 - j].gen;

            MPI_Sendrecv(send, pair, MPI_UNSIGNED, to, tag1, recv, pair, MPI_UNSIGNED,
                from, tag2, MPI_COMM_WORLD, &status);

            for (uint j = 0; j < pair; j++)
                generat[j].gen = recv[j];

            delete [] send;
            delete [] recv;
        }
    }

    for (uint j = 0; j < GEN_BEGIN; j++)
        generat[j].gen_fit = fitness(generat[j], items, volume);
    sort(generat.begin(), generat.end(), comp);

    return;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
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
        ifstream infile(argv[1]);
        if (infile.is_open() == false)
        {
            cerr << ">Can not open infile with such name." << endl;
            return -1;
        }

        srand(time(0));

        int volume = 0;
        int quan = 0;
        vector< item > items;

        infile >> volume;
        infile >> quan;

        items.resize(quan);

        for (uint i = 0; i < quan; i++)
        {
            infile >> items[i].vol; /*volume*/
            infile >> items[i].val; /*value*/
        }

        for (uint i = 1; i < nProc; i++)
        {
            int rand_val;
            rand_val = rand();
            MPI_Send(&rand_val, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
        }
        
        MPI_Bcast(&volume, 1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(&quan, 1, MPI_INT, root, MPI_COMM_WORLD);
        int *tmp_vol = new int [quan];
        int *tmp_val = new int [quan];

        for (uint i = 0; i < quan; i++)
        {
            tmp_vol[i] = items[i].vol;
            tmp_val[i] = items[i].val;
        }

        MPI_Bcast(tmp_vol, quan, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(tmp_val, quan, MPI_INT, root, MPI_COMM_WORLD);

        delete [] tmp_vol;
        delete [] tmp_val;

        vector<generation> generat;

        generat.resize(GEN_BEGIN);
        int top = 1;
        top = (top << quan) - 1;

        for (uint i = 0; i < GEN_BEGIN; i++)
            generat[i].gen = (int)(rand()) % top + 1;

        selection(items, generat, volume, 1, nProc - 1);

        std::vector<uint> res_vec;
        res_vec.resize(nProc);
        res_vec[0] = generat[GEN_BEGIN - 1].gen;

        for (uint i = 1; i < nProc; i++)
        {
            uint gen;
            MPI_Recv(&gen, 1, MPI_UNSIGNED, i, tag, MPI_COMM_WORLD, &status);
            res_vec[i] = gen;
        }

        sort(res_vec.begin(), res_vec.end());

        uint one = 1 << (quan - 1), count = quan - 1;
        //uint tmp = 0, fit = 0, v = 0, res = generat[GEN_BEGIN - 1].gen;
        uint tmp = 0, fit = 0, v = 0, res = res_vec[nProc - 1];

        cout << endl << ">Best chrom:" << endl;
        for (int i = quan - 1; i > -1; i--)
        {
            tmp = res & one;
            tmp = tmp >> count;
            fit +=  tmp * items[quan - 1 - i].val;
            v += tmp * items[quan - 1 - i].vol;
            cout << tmp << " ";
            count--;
            one = one >> 1;
        }
        cout << endl << endl << ">Value = " << fit << endl << ">Volume = " << v << endl;

        infile.close();
    }
    else
    {
        int rand_val = 0;
        MPI_Recv(&rand_val, 1, MPI_INT, root, tag, MPI_COMM_WORLD, &status);

        srand(rand_val);

        int quan, volume;
        MPI_Bcast(&volume, 1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(&quan, 1, MPI_INT, root, MPI_COMM_WORLD);

        int *tmp_vol = new int [quan];
        int *tmp_val = new int [quan];
        MPI_Bcast(tmp_vol, quan, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(tmp_val, quan, MPI_INT, root, MPI_COMM_WORLD);

        vector< item > items;
        items.resize(quan);

        for (uint i = 0; i < quan; i++)
        {
            items[i].vol = tmp_vol[i]; /*volume*/
            items[i].val = tmp_val[i]; /*value*/
        }

        delete [] tmp_vol;
        delete [] tmp_val;

        vector<generation> generat;

        generat.resize(GEN_BEGIN);
        int top = 1;
        top = (top << quan) - 1;

        for (uint i = 0; i < GEN_BEGIN; i++)
            generat[i].gen = (int)(rand()) % top + 1;

        selection(items, generat, volume, (myRank + 1) % nProc, myRank - 1);

        MPI_Send(&generat[GEN_BEGIN - 1].gen, 1, MPI_UNSIGNED, root, tag, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}