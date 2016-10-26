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

int main(int argc, char const *argv[])
{
    if (argc != 3) /*make it !!4!! after debug*/
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    } 

    /*ifstream outfile(argv[3]);
    if (outfile.is_open() == false)
    {
        cerr << ">Can not open outfile with such name." << endl;
        return -1;
    }*/

    uint left = atoi(argv[1]), right = atoi(argv[2]);

    if (right <= left)
    {
        cerr << ">Wrong arguments: left >= right." << endl;
        return -1;
    }

    std::vector<bool> primes(right - left + 2, 1);
    primes[0] = primes[1] = 0;

    for (uint i = left + 1; i * i < right + 1; i++)
        if (primes[i] == 1)
            for (int j = i * i; j < right + 1; j += i)
                primes[j] = 0;

    for (uint i = 0; i < primes.size(); i++)
        if (primes[i] == 1)
            cout << i << endl;

    //outfile.close();
    return 0;
}