#include <iostream>
#include <ctime> // в ней функция time
#include <cstdlib>
#include <fstream>
 
using namespace std;
 
int main()
{
    ofstream matrix("matrixB_1024x4096");
    srand(time(NULL)); // Инициализируем генератор случайных чисел. 
    int n = 0, m = 0; 
    cin >> n; // Считываем с клавиатуры n строки
    cin >> m; // Считываем с клавиатуры m столбцы
    matrix << n << " " << m << endl;
    double **a = new double* [n]; // Создаем массив указателей
    for (int i = 0; i < n; i++)
    {
        a[i] = new double [m]; // Создаем элементы
    }
    // А дальше работа как с обычным массивом. 
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            //a[i][j] = rand() % 100; // Каждый элемент случайному числу от 0 до 9
            a[i][j] =  (double)(rand()) / RAND_MAX * 200 - 100;
            matrix << a[i][j] << " "; // Вывести элементы на консольку
        }
        matrix << endl; // Двумерный массив. Строка кончилась, переводим строку и на консоли
    }
    // Удаление массива
    for (int i = 0; i < n; i++)
    {
        delete[]a[i]; // Удаляем каждый элемент
    }
    delete [] a; // А потом массив
    matrix.close();
    return 0;
}