#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include "mpi.h"
#define EPS 0.0001
using namespace std;
int main(int argc, char **argv)
{

	MPI_Init(&argc, &argv);
	int proc = 0, myRank, n;
	MPI_Comm_size(MPI_COMM_WORLD, &proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	cout.precision(10);
	MPI_File fh_a;
	MPI_File fh_b;
	MPI_File fh_c;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_a);
	MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_b);
	MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_c);
	
	double *A, *B, *C;
	n = atoi(argv[4]);

	A = (double*)malloc(sizeof(double) * n * n);
	B = (double*)malloc(sizeof(double) * n * n);
	C = (double*)malloc(sizeof(double) * n * n);

	MPI_File_seek(fh_a, sizeof(int), MPI_SEEK_SET);
	MPI_File_seek(fh_b, sizeof(int), MPI_SEEK_SET);
	
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
		{
			double tmp;
			MPI_File_read(fh_a, &tmp, 1, MPI_DOUBLE, &status);
			A[i*n + j] = tmp;
			MPI_File_read(fh_b, &tmp, 1, MPI_DOUBLE, &status);	
			B[i*n + j] = tmp;
			//cout << tmp << endl;
		}

	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
		{
			double res = 0.0;
			for(int k = 0; k < n; k++)
			{
				res += A[i*n + k] * B[k*n + j];
			}
			C[i*n + j] = res;
			//cout << res << endl;
		}

	MPI_File_seek(fh_c, sizeof(int), MPI_SEEK_SET);
	for(int i = 0; i < n*n; i++)
	{
		double tmp;
		MPI_File_read(fh_c, &tmp, 1, MPI_DOUBLE, &status);
		//cout << tmp << "   " << C[i] << endl;
		if(fabs(fabs(tmp) - fabs(C[i])) > EPS)
		{
			cout << "ERRROR " << tmp << "  " << C[i] << endl;
			return 0;
		}

	}
	cout << "OK!" << endl;
	MPI_File_close(&fh_a);
	MPI_File_close(&fh_b);
	MPI_File_close(&fh_c);
	MPI_Finalize();
	return 0;
}
	