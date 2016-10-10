#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

#include <vector>
#include <algorithm>

#define GEN_BEGIN 300
#define GEN_END 1000
#define PERS 0.3
#define MUT_PERS 3

using namespace std;

struct item
{
    double vol;
    double val;
};

struct parents
{
    uint par1;
    uint par2;
};

struct generation
{
    uint gen;
    double gen_fit;
};

struct scomp
{
    bool operator () (generation fit1, generation fit2)
    {
        return fit1.gen_fit < fit2.gen_fit;
    }
} comp;

void mutation (uint &gen, uint &top)
{
    uint pos = 0, point = 1;

    pos = (uint)(rand()) % top + 1;

    point = point << pos;

    gen = gen | point; /*???????????*/

    return;
}

int fitness (const generation &individual, const vector< item > &items, const int &volume)
{
    int fit = 0, v = 0;
    uint one = 1, tmp = 0, conut = 0;
    uint quan = items.size();

    for (int i = quan - 1; i > -1; i--)
    {
        tmp = individual.gen & one;
        tmp = tmp >> conut;
        fit +=  tmp * items[i].val;
        v += tmp * items[i].vol;
        one = one << 1;
        conut++;
    }

    if (v > volume)
        fit = 0;

    if (fit != 0)
        cout << fit << " " << v << endl;

    return fit;
}

uint cross (parents par, uint &quan)
{
    uint point = 0, res = 0;
    uint top = 1;
    top = (top << quan) - 1;

    point = (uint)(rand()) % top + 1;

    par.par1 = par.par1 >> point;
    par.par2 = par.par2 << point;

    par.par1 = par.par1 << point;
    par.par2 = par.par2 >> point;

    res = par.par1 | par.par2;

    point = (uint)(rand()) % 100 + 1;
    if (point < MUT_PERS)
        mutation(res, top);

    return res;
}

void selection (const vector< item > &items, vector<generation> &generat, const int &volume)
{
    vector<parents> parent;
    uint quan = items.size();

    uint pair = GEN_BEGIN * PERS;

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
    }

    for (uint j = 0; j < GEN_BEGIN; j++)
        generat[j].gen_fit = fitness(generat[j], items, volume);
    sort(generat.begin(), generat.end(), comp);

    return;
}

int main(int argc, char const *argv[])
{
    if (argc != 2)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    } 

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

    vector<generation> generat;

    generat.resize(GEN_BEGIN);
    int top = 1;
    top = (top << quan) - 1;

    for (uint i = 0; i < GEN_BEGIN; i++)
        generat[i].gen = (int)(rand()) % top + 1;

    selection(items, generat, volume);

    cout << generat[GEN_BEGIN - 1].gen << endl;

    fitness(generat[GEN_BEGIN - 1], items, volume);

    infile.close();
    return 0;
}