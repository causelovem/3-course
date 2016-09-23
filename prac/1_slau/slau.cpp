#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include <vector>
#include <algorithm>

#define EPS 0.000000000001

using namespace std;

bool flag = false;

void make_solution (std::vector < std::vector <double> > &matrix, std::vector <double> &right_vec, std::vector <int> var_vec, int num_eq, int num_var);
std::vector <int> new_gauss (std::vector < std::vector <double> > &matrix, std::vector <double> &right_vec, int num_eq, int num_var);
std::vector <double> gauss (std::vector < std::vector <double> > matrix, std::vector <double> right_vec, int num_eq, int num_var);

int main(int argc, char const *argv[])
{
    double x, y;
    x = (float)clock();

    std::vector < std::vector <double> > matrix, tmp_matrix;
    std::vector <double> tmp_vec, right_vec;
    std::vector <int> var_vec;

    int number_of_equesions = 0, number_of_variable = 0;
    double num = 0;

    cout << ">Input number of equesions\n";
    cin >> number_of_equesions;

    cout << ">Input number of variable\n";
    cin >> number_of_variable;

    cout << ">Input matrix\n";

    for (int i = 0; i < number_of_equesions; i++)
    {
        for (int j = 0; j < number_of_variable; j++)
        {
            cin >> num;
            tmp_vec.push_back(num);
        }
        matrix.push_back(tmp_vec);
        tmp_vec.clear();
    }

    cout << ">Input right parts\n";

    for (int i = 0; i < number_of_equesions; i++)
    {
        cin >> num;
        right_vec.push_back(num);
    }

    cout << ">Your slau\n";
    for (int i = 0; i < number_of_equesions; i++)
    {
        for (int j = 0; j <  number_of_variable; j++)
        {
            cout << matrix[i][j] << "*x" << j;
            if (j != number_of_variable - 1)
                cout << " + ";
        }
        cout << " = " << right_vec[i];
        cout << endl;
    }

    tmp_vec = right_vec;
    tmp_matrix = matrix;

    var_vec = new_gauss(tmp_matrix, tmp_vec, number_of_equesions, number_of_variable);

    /*for (int t = 0; t < number_of_equesions; t++)
    {
        for (int j = 0; j < number_of_variable; j++)
            cout << tmp_matrix[t][j] << " ";

        cout << " " << tmp_vec[t];
        cout << endl;
    }*/

    if (flag == true)
    {
        cout << ">There is no solution for your slau :(" << endl;

        y = (float)clock();
        cout << endl << ">TIME " << (float)(y - x) / CLOCKS_PER_SEC << endl;

        return -1;
    }

    make_solution(tmp_matrix, tmp_vec, var_vec, number_of_equesions, number_of_variable);

    y = (float)clock();
    cout << endl << ">TIME " << (float)(y - x) / CLOCKS_PER_SEC << endl;

    return 0;
}

std::vector <double> gauss (std::vector < std::vector <double> > matrix, std::vector <double> right_vec, int num_eq, int num_var)
{
    std::vector <double> tmp_vec;

    int step = 0, max_num = 0;
    double max = 0, tmp = 0;

    while (step < num_eq)
    {
        tmp_vec.clear();
        max = fabs(matrix[step][step]);
        max_num = step;

        for (int i = step + 1; i < num_eq; i++)
        {
            if (fabs(matrix[i][step]) > max)
            {
                max = fabs(matrix[i][step]);
                max_num = i;
            }
        }

        if (((max < EPS) && (step != num_eq - 1)) || ((max < EPS) && (step == num_eq - 1) && right_vec[step]))
        /*if (max < EPS)*/
        {
            cout << ">Can not solve slau because of zero col/row\n";
            return tmp_vec;
        }
        else
        if ((max < EPS) && (step == num_eq - 1))
            flag = true;

        tmp_vec = matrix[step];
        matrix[step] = matrix[max_num];
        matrix[max_num] = tmp_vec;

        tmp = right_vec[step];
        right_vec[step] = right_vec[max_num];
        right_vec[max_num] = tmp;

        for (int i = step; i < num_eq; i++)
        {
            tmp = matrix[i][step];
            if (fabs(tmp) < EPS)
                continue;
            for (int j = 0; j < num_eq; j++)
                matrix[i][j] = matrix[i][j] / tmp;

            right_vec[i] = right_vec[i] / tmp;

            if (i == step)
                continue;

            for (int j = 0; j < num_eq; j++)
                matrix[i][j] = matrix[i][j] - matrix[step][j];


            right_vec[i] = right_vec[i] - right_vec[step];
        }

        /*for (int i = 0; i < num_eq; i++)
        {
            for (int j = 0; j < num_var; j++)
                cout << matrix[i][j] << " ";

            cout << " " << right_vec[i];
            cout << endl;
        }*/

        step++;
    }

    tmp_vec.clear();
    for (int i = num_var - 1; i > -1; i--)
    {
        tmp_vec.push_back(right_vec[i]);
        for (int j = 0; j < i; j++)
            right_vec[j] = right_vec[j] - matrix[j][i] * right_vec[i];
    }

    return tmp_vec;
}

std::vector <int> new_gauss (std::vector < std::vector <double> > &matrix, std::vector <double> &right_vec, int num_eq, int num_var)
{
    std::vector <double> tmp_vec;
    std::vector <int> var_vec;

    int lim = 0, col = 0, row = 0, i = 0;
    double tmp = 0, max_num = matrix[0][0];
    bool break_flag = false;

    if (num_eq >= num_var)
        lim = num_var;
    else
        lim = num_eq;

    for (i = 0; i < num_var; i++)
        var_vec.push_back(i);

    for (i = 0; i < lim; i++)
    {
        tmp_vec.clear();
        max_num = matrix[i][i];
        row = i;
        col = i;
        for (int k = i; k < num_var; k++)
            for (int l = i; l < num_eq; l++)
                if (fabs(matrix[l][k]) > fabs(max_num))
                {
                    max_num = matrix[l][k];
                    row = l;
                    col = k;
                }

        tmp_vec = matrix[row];
        matrix[row] = matrix[i];
        matrix[i] = tmp_vec;

        for (int j = 0; j < num_eq; j++)
        {
            tmp = matrix[j][col];
            matrix[j][col] = matrix[j][i];
            matrix[j][i] = tmp;
        }

        tmp = right_vec[row];
        right_vec[row] = right_vec[i];
        right_vec[i] = tmp;

        tmp = var_vec[col];
        var_vec[col] = var_vec[i];
        var_vec[i] = tmp;

        if (fabs(max_num) < EPS)
        {
            break_flag = true;
            break;
        }

        /*for (int t = 0; t < num_eq; t++)
        {
            for (int j = 0; j < num_var; j++)
                cout << matrix[t][j] << " ";

            cout << " " << right_vec[t];
            cout << endl;
        }*/

        for (int k = i + 1; k < num_eq; k++)
            if (fabs(matrix[k][i]) > EPS)
            {
                tmp = matrix[k][i];
                //cout << "&";
                for (int l = i; l < num_var; l++)
                {
                    matrix[k][l] = matrix[k][l] * max_num / tmp;
                    if (fabs(matrix[k][l]) < EPS)
                        matrix[k][l] = 0;
                }

                for (int l = i; l < num_var; l++)
                {
                    matrix[k][l] = matrix[k][l] - matrix[i][l];
                    if (fabs(matrix[k][l]) < EPS)
                        matrix[k][l] = 0;
                }
                
                right_vec[k] = right_vec[k] * max_num / tmp;
                right_vec[k] = right_vec[k] - right_vec[i];
                if (fabs(right_vec[k]) < EPS)
                    right_vec[k] = 0;
            }

        /*for (int t = 0; t < num_eq; t++)
        {
            for (int j = 0; j < num_var; j++)
                cout << matrix[t][j] << " ";

            cout << " " << right_vec[t];
            cout << endl;
        }*/
    }

    tmp_vec.clear();

    if ((break_flag == true) || (lim == num_var))
        for (int j = i; j < num_eq; j++)
            if (right_vec[j] != 0)
            {
                flag = true;
                return var_vec;
            }

    for (int k = 0; k < i; k++)
    {
        if (fabs(matrix[k][k]) > EPS)
        {
            tmp = matrix[k][k];
            for (int l = k; l < num_var; l++)
            {
                matrix[k][l] = matrix[k][l] / tmp;
                if(fabs(matrix[k][l]) < EPS)
                    matrix[k][l] = 0;
            }

            right_vec[k] = right_vec[k] / tmp;
            if (fabs(right_vec[k]) < EPS)
                right_vec[k] = 0;
        }

        /*for (int t = 0; t < num_eq; t++)
        {
            for (int j = 0; j < num_var; j++)
                cout << matrix[t][j] << " ";

            cout << " " << right_vec[t];
            cout << endl;
        }*/
    }

    return var_vec;
}

void make_solution (std::vector < std::vector <double> > &matrix, std::vector <double> &right_vec, std::vector <int> var_vec, int num_eq, int num_var)
{
    int lim = 0, step = 0, tmp_int = 0;
    double tmp_doub = 0;
    std::vector <double> tmp_vec;

    if (num_eq >= num_var)
        lim = num_var;
    else
        lim = num_eq;

    for (step = 0; step < lim; step++)
        if (fabs(matrix[step][step]) == 0)
        {
            break;
        }
        step--;

    //cout << step << endl;

    if ((step == lim - 1) && (lim == num_var))
    {
        for (int k = step; k > -1; k--)
            for (int l = 0; l < k; l++)
            {
                right_vec[l] = right_vec[l] - matrix[l][k] * right_vec[k];
                if (fabs(right_vec[l]) < EPS)
                    right_vec[l] = 0;
            }

        for (int k = 0; k < num_var - 1; k++)
            for (int l = k + 1; l < num_var; l++)
                if (var_vec[k] > var_vec[l])
                {
                    tmp_int = var_vec[k];
                    var_vec[k] = var_vec[l];
                    var_vec[l] = tmp_int;

                    tmp_doub = right_vec[k];
                    right_vec[k] = right_vec[l];
                    right_vec[l] = tmp_doub;
                }

        cout << "Your solution *_*" << endl;
        for (int k = 0; k < num_var; k++)
            cout << "x" << k + 1 << " = " << right_vec[k] << "   ";
        cout << endl;
    }
    else
    {
        cout << ">Input next variable:" << endl;
        for (int k = step + 1; k < num_var; k++)
        {
            cout << "x" << var_vec[k] + 1 << " = ";
            cin >> right_vec[k];
        }

        for (int k = num_var - 1; k > -1; k--)
            for (int l = 0; l < k; l++)
            {
                right_vec[l] = right_vec[l] - matrix[l][k] * right_vec[k];
                if (fabs(right_vec[l]) < EPS)
                    right_vec[l] = 0;
                
                /*for (int t = 0; t < num_eq; t++)
                {
                    for (int j = 0; j < num_var; j++)
                        cout << matrix[t][j] << " ";

                    cout << " " << right_vec[t];
                    cout << endl;
                }*/
            }

        for (int k = 0; k < num_var - 1; k++)
            for (int l = k + 1; l < num_var; l++)
                if (var_vec[k] > var_vec[l])
                {
                    tmp_int = var_vec[k];
                    var_vec[k] = var_vec[l];
                    var_vec[l] = tmp_int;

                    tmp_doub = right_vec[k];
                    right_vec[k] = right_vec[l];
                    right_vec[l] = tmp_doub;
                }

        cout << "Your solution *_*" << endl;
        for (int k = 0; k < num_var; k++)
            cout << "x" << k + 1 << " = " << right_vec[k] << "   ";
        cout << endl;
    }

    return;
}