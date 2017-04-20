#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <omp.h>
#include <complex>

using namespace std;

#define EPS 0.001
 
typedef std::complex <double> complexd;

int main(int argc, char const *argv[])
{
    if (argc != 5)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    ifstream vecAfile(argv[1]);
    if (vecAfile.is_open() == false)
    {
        cerr << ">Can not open vecA with such name." << endl;
        return -1;
    }

    ifstream vecBfile(argv[2]);
    if (vecBfile.is_open() == false)
    {
        cerr << ">Can not open vecB with such name." << endl;
        vecAfile.close();
        return -1;
    }

    int n = 0, k = atoi(argv[3]), size = 0;

    vecAfile >> n;
    vecBfile >> n;

    if (k > n)
    {
        cerr << ">Number of bit is too big (K > N)." << endl;
        vecAfile.close();
        vecBfile.close();
        return -1;
    }

    omp_set_num_threads(atoi(argv[4]));

    size = 1 << n;     /*2 ^ n*/

    complexd *vecA = NULL, *vecB = NULL;

    vecA = new complexd [size];
    vecB = new complexd [size];

    double re = 0, im = 0;

    for (int i = 0; i < size; i++)
    {
        vecAfile >> re;
        vecAfile >> im;

        complexd tmp1 (re, im);
        vecA[i] = tmp1;

        vecBfile >> re;
        vecBfile >> im;

        complexd tmp2 (re, im);
        vecB[i] = tmp2;
    }

    /*
                      (1  1)
    H = (1 / sqrt(2)) (1 -1)

    */

    int step = 1 << (n - k + 1), loop1 = 1 << (k - 1), loop2 = 1 << (n - k), mask = 1 << (n - k);

    //float time = (float)clock();

    /*#pragma omp parallel for
    for (int j = 0; j < loop1; j++)
    {
        int tmp = j * step;

        #pragma omp parallel for
        for (int i = 0; i < loop2;  i++)
        {
            complexd tmp1 = (1 / sqrt(2)) * (vecB[i + tmp] + vecB[i + tmp + loop2]);
            complexd tmp2 = (1 / sqrt(2)) * (vecB[i + tmp] - vecB[i + tmp + loop2]);
            vecB[i + tmp] = tmp1;
            vecB[i + tmp + loop2] = tmp2;
        }
    }*/

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        int ormask = (i | mask);
        complexd tmp1 = (1 / sqrt(2)) * (vecB[i] + vecB[ormask]);
        complexd tmp2 = (1 / sqrt(2)) * (vecB[i] - vecB[ormask]);

        //if (ormask > andmask)
        if (i < ormask)
        {
            vecB[i] = tmp1;
            vecB[ormask] = tmp2;
        }
    }

    //cout << ">Time of computation = " << ((float)clock() - time) / CLOCKS_PER_SEC << endl;
    cout << "N = " << n << "\nK = " << k << endl;
    cout << ">Number of threads = " << atoi(argv[4]) << endl;

    for (int i = 0; i < size; i++)
        if ((abs(vecA[i]) - abs(vecB[i])) > EPS)
        {
            cout << "NO :(" << endl;
            delete [] vecA;
            delete [] vecB;
            vecAfile.close();
            vecBfile.close();

            return -1;
        }
    cout << "YES *_*\nvecA == vecB" << endl;

    vecAfile.close();
    vecBfile.close();

    delete [] vecA;
    delete [] vecB;

    return 0;
}