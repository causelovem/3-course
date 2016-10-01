#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void yakobi (double **matrix, int size, ofstream &out_file);
void mul_matrix (double **&matrixA, double **&matrixB, int size);

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        printf(">Unexpected quantity of arguments, check your comand string.\n");
        return -1;
    }

    ofstream out_file(argv[2]);
    ifstream matrix_file(argv[1]);

    if (matrix_file.is_open() == false)
    {
        printf(">Can not open matrix file with such name.\n");
        return -1;
    }

    double **matrix =  NULL;
    int size = 0;

    matrix_file >> size;

    matrix = new double* [size];
    for (int i = 0; i < size; i++)
        matrix[i] = new double [size];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[i][j] = 0;

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix_file >> matrix[i][j];

    yakobi(matrix, size, out_file);

    for (int i = 0; i < size; i++)
    {
        delete [] matrix[i];
        matrix[i] = NULL;
    }
    delete [] matrix;
    matrix = NULL;

    matrix_file.close();
    out_file.close();
    return 0;
}

void yakobi (double **matrix, int size, ofstream &out_file)
{
    double eps = 0.001;
    double **vec_matrix = NULL;
    double **rot_matrix = NULL;
    double max_num = 0, tmp = 0, rsin = 0, rcos = 0;
    int row = 0, col = 0;
    bool flag = false;

    vec_matrix = new double* [size];
    for (int i = 0; i < size; i++)
        vec_matrix[i] = new double [size];

    rot_matrix = new double* [size];
    for (int i = 0; i < size; i++)
        rot_matrix[i] = new double [size];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            vec_matrix[i][j] = 0;

    while(1)
    {
        max_num = matrix[0][1];
        row = 0;
        col = 1;
        for (int i = 1; i < size; i++)
            for (int j = 0; j < i; j++)
                if (fabs(matrix[i][j]) > fabs(max_num))
                {
                    max_num = matrix[i][j];
                    row = i;
                    col = j;
                }

        if (fabs(max_num) < eps)
            break;

        tmp = 2 * matrix[row][col] / (matrix[row][row] - matrix[col][col]);

        //cout << (tmp / fabs(tmp)) << endl;

        rsin = sqrt((0.5) * (1 - 1 / (sqrt(1 + tmp * tmp))));
        rcos = sqrt((0.5) * (1 + 1 / (sqrt(1 + tmp * tmp))));

        //cout << tmp << " " << rsin << " " << rcos << endl;

        //break;

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                rot_matrix[i][j] = 0;
            rot_matrix[i][i] = 1;
        }

        rot_matrix[row][row] = rcos;
        rot_matrix[row][col] = -rsin;
        rot_matrix[col][row] = rsin;
        rot_matrix[col][col] = rcos;

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                cout << rot_matrix[j][i] << " ";
            }

            cout << endl;
        }

        if (flag == true)
            mul_matrix(vec_matrix, rot_matrix, size);
        else
        {
            flag = true;

            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    vec_matrix[i][j] = rot_matrix[i][j];
        }

        mul_matrix(matrix, rot_matrix, size);

        rot_matrix[row][row] = rcos;
        rot_matrix[row][col] = rsin;
        rot_matrix[col][row] = -rsin;
        rot_matrix[col][col] = rcos;

        mul_matrix(rot_matrix, matrix, size);

        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                matrix[i][j] = rot_matrix[i][j];

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                cout << matrix[j][i] << " ";
            }

            cout << endl;
        }

        break;
    }

    /*for (int i = 0; i < size; i++)
        cout << matrix[i][i] << " ";
    cout << endl;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << vec_matrix[j][i] << " ";
        }

        cout << endl;
    }*/

    for (int i = 0; i < size; i++)
        out_file << matrix[i][i] << " ";
    out_file << endl;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            out_file << vec_matrix[j][i] << " ";
        }

        out_file << endl;
    }

    for (int i = 0; i < size; i++)
    {
        delete [] vec_matrix[i];
        vec_matrix[i] = NULL;
        delete [] rot_matrix[i];
        rot_matrix[i] = NULL;
    }
    delete [] vec_matrix;
    vec_matrix = NULL;
    delete [] rot_matrix;
    rot_matrix = NULL;

    return;
}

void mul_matrix (double **&matrixA, double **&matrixB, int size)
{
    double **matrixC = NULL;

    matrixC = new double* [size];
    for (int i = 0; i < size; i++)
        matrixC[i] = new double [size];

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrixC[i][j] = 0; 

    for (int i = 0; i < size; i++)
    {
        for (int k = 0; k < size; k++)
        {
            double elem = matrixA[i][k];
            for (int j = 0; j < size; j++)
                matrixC[i][j] += elem * matrixB[k][j];
        }
    }

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrixA[i][j] = matrixC[i][j];

    for (int i = 0; i < size; i++)
    {
        delete [] matrixC[i];
        matrixC[i] = NULL;
    }
    delete [] matrixC;
    matrixC = NULL;

    return;
}