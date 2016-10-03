#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void yakobi (double **&matrix, int size, ofstream &out_file);
void mul_matrix (double **&matrixA, double **&rot_matrix, int size, int i_cord, int j_cord, bool flag);

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
        out_file.close();
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

void yakobi (double **&matrix, int size, ofstream &out_file)
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

        rsin = (tmp / fabs(tmp)) * sqrt((0.5) * (1.0 - 1.0 / (sqrt(1.0 + tmp * tmp))));
        rcos = sqrt((0.5) * (1.0 + 1.0 / (sqrt(1.0 + tmp * tmp))));

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

        if (flag == true)
            mul_matrix(vec_matrix, rot_matrix, size, row, col, true);
        else
        {
            flag = true;

            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    vec_matrix[i][j] = rot_matrix[i][j];
        }

        mul_matrix(matrix, rot_matrix, size, row, col, false);
    }

    for (int i = 0; i < size; i++)
        out_file << "l" << i + 1 << " = " << matrix[i][i] << ";  ";
    out_file << endl;

    for (int i = 0; i < size; i++)
    {
        out_file << "v" << i + 1 << " = { "; 
        for (int j = 0; j < size; j++)
        {
            out_file << vec_matrix[j][i] << "; ";
        }

        out_file << "}"<< endl;
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

void mul_matrix (double **&matrixA, double **&rot_matrix, int size, int i_cord, int j_cord, bool flag)
{
    double *tmp_vec1 = NULL, *tmp_vec2 = NULL;

    tmp_vec1 = new double [size];
    tmp_vec2 = new double [size];

    for (int i = 0; i < size; i++)
    {
        tmp_vec1[i] = 0;
        tmp_vec2[i] = 0;
    }

    for (int i = 0; i < size; i++)
    {
        tmp_vec1[i] = matrixA[i][i_cord] * rot_matrix[i_cord][i_cord] + matrixA[i][j_cord] * rot_matrix[j_cord][i_cord];
        tmp_vec2[i] = matrixA[i][i_cord] * rot_matrix[i_cord][j_cord] + matrixA[i][j_cord] * rot_matrix[j_cord][j_cord];
    }

    for (int i = 0; i < size; i++)
    {
        matrixA[i][i_cord] = tmp_vec1[i];
        matrixA[i][j_cord] = tmp_vec2[i];
    }

    for (int i = 0; i < size; i++)
    {
        tmp_vec1[i] = 0;
        tmp_vec2[i] = 0;
    }

    if (flag == true)
    {
        delete [] tmp_vec1;
        delete [] tmp_vec2;
        tmp_vec1 = NULL;
        tmp_vec2 = NULL;
        return;
    }

    rot_matrix[i_cord][j_cord] *= -1;
    rot_matrix[j_cord][i_cord] *= -1;

    for (int i = 0; i < size; i++)
    {
        tmp_vec1[i] = matrixA[i_cord][i] * rot_matrix[i_cord][i_cord] + matrixA[j_cord][i] * rot_matrix[i_cord][j_cord];
        tmp_vec2[i] = matrixA[i_cord][i] * rot_matrix[j_cord][i_cord] + matrixA[j_cord][i] * rot_matrix[j_cord][j_cord];
    }

    for (int i = 0; i < size; i++)
    {
        matrixA[i_cord][i] = tmp_vec1[i];
        matrixA[j_cord][i] = tmp_vec2[i];
    }

    delete [] tmp_vec1;
    delete [] tmp_vec2;
    tmp_vec1 = NULL;
    tmp_vec2 = NULL;

    return;
}