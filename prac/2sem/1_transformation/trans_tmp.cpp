#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <omp.h>
#include <complex>

using namespace std;

#define EPS 0.000001
 
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

    ofstream vecBfile(argv[2]);
    if (vecBfile.is_open() == false)
    {
        cerr << ">Can not open vecB with such name." << endl;
        vecAfile.close();
        return -1;
    }

    int n = 0, k = atoi(argv[3]), size = 0;

    vecAfile >> n;

    if (k > n)
    {
        cerr << ">Number of bit is too big (K > N)." << endl;
        vecAfile.close();
        vecBfile.close();
        return -1;
    }

    omp_set_num_threads(atoi(argv[4]));

    size = 1 << n;     /*2 ^ n*/

    complexd *vecA = NULL, *vecB = NULL, *vecC = NULL;

    vecA = new complexd [size];
    vecB = new complexd [size];
    vecC = new complexd [size];

    double re = 0, im = 0;

    for (int i = 0; i < size; i++)
    {
        vecAfile >> re;
        vecAfile >> im;

        complexd tmp (re, im);
        vecA[i] = tmp;
        vecC[i] = tmp;
    }

    /*cout << "VECTOR A" << endl;
    for (int i = 0; i < size; i++)
        cout << vecA[i] << endl;
    cout << endl;*/

    /*
                      (1  1)
    H = (1 / sqrt(2)) (1 -1)

    */

    int step = 1 << (n - k + 1), loop1 = 1 << (k - 1), loop2 = 1 << (n - k);

    float time = (float)clock();

    #pragma omp parallel for
    for (int j = 0; j < loop1; j++)
    {
        int tmp = j * step;

        #pragma omp parallel for
        for (int i = 0; i < loop2;  i++)
        {
            /*vecB[i + tmp] = (1 / sqrt(2)) * (vecA[i + tmp] + vecA[i + tmp + loop2]);
            vecB[i + tmp + loop2] = (1 / sqrt(2)) * (vecA[i + tmp] - vecA[i + tmp + loop2]);*/

            complexd tmp1 = (1 / sqrt(2)) * (vecA[i + tmp] + vecA[i + tmp + loop2]);
            complexd tmp2 = (1 / sqrt(2)) * (vecA[i + tmp] - vecA[i + tmp + loop2]);
            vecA[i + tmp] = tmp1;
            vecA[i + tmp + loop2] = tmp2;
        }
    }

    cout << ">Time of computation = " << ((float)clock() - time) / CLOCKS_PER_SEC << endl;
    cout << "N = " << n << "\nK = " << k << endl;
    cout << ">Number of threads = " << atoi(argv[4]) << endl;

    /*cout << "VECTOR B" << endl;
    for (int i = 0; i < size; i++)
        cout << vecB[i] << endl;*/

    #pragma omp parallel for
    for (int j = 0; j < loop1; j++)
    {
        int tmp = j * step;

        #pragma omp parallel for
        for (int i = 0; i < loop2;  i++)
        {
            /*vecC[i + tmp] = (1 / sqrt(2)) * (vecB[i + tmp] + vecB[i + tmp + loop2]);
            vecC[i + tmp + loop2] = (1 / sqrt(2)) * (vecB[i + tmp] - vecB[i + tmp + loop2]);*/

            complexd tmp1 = (1 / sqrt(2)) * (vecA[i + tmp] + vecA[i + tmp + loop2]);
            complexd tmp2 = (1 / sqrt(2)) * (vecA[i + tmp] - vecA[i + tmp + loop2]);
            vecA[i + tmp] = tmp1;
            vecA[i + tmp + loop2] = tmp2;
        }
    }

    /*cout << endl <<  "VECTOR C (just for check)" << endl;
    for (int i = 0; i < size; i++)
        cout << vecC[i] << endl;*/

    for (int i = 0; i < size; i++)
        if ((abs(vecA[i]) - abs(vecC[i])) > EPS)
        {
            cout << "NO :(" << endl;
            delete [] vecA;
            delete [] vecB;
            delete [] vecC;
            vecAfile.close();
            vecBfile.close();

            return -1;
        }
    cout << "YES *_*\nvecA == vecC" << endl;

    /*vecBfile << n << endl;
    for (int i = 0; i < size; i++)
    {
        vecBfile << double(real(vecB[i])) << " ";
        vecBfile << double(imag(vecB[i])) << endl;
    }*/

    vecAfile.close();
    vecBfile.close();

    delete [] vecA;
    delete [] vecB;
    delete [] vecC;

    return 0;
}