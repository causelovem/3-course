#include <cstdlib>
#include <iostream>
#include <fstream>

#include <omp.h>

using namespace std;

int main(int argc, char const *argv[])
{
	if (argc != 4)
    {
        cerr << ">Unexpected quantity of arguments, check your comand string." << endl;
        return -1;
    }

    ofstream vecAfile(argv[1]);
    if (vecAfile.is_open() == false)
    {
        cerr << ">Can not open vecA with such name." << endl;
        return -1;
    }

    int n = atoi(argv[2]), size = 0;

    omp_set_num_threads(atoi(argv[3]));

    srand(time(0));

    size = 1 << n;     /*2 ^ n*/

    vecAfile << n << endl;

    //#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
    	vecAfile << (double)(rand()) / RAND_MAX * 20 - 10 << " ";
    	vecAfile << (double)(rand()) / RAND_MAX * 20 - 10 << endl;
    }

    vecAfile.close();
	return 0;
}