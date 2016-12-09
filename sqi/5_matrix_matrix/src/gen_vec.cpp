#include <iostream>
#include <ctime> // в ней функция time
#include <cstdlib>
#include <fstream>
 
using namespace std;
 
int main()
{
    ofstream matrix("matrixB_4096");
    srand(time(NULL)); // Инициализируем генератор случайных чисел. 
    int n = 0; 
    cin >> n; // Считываем с клавиатуры n
    matrix << n << endl;
    double *a = new double [n]; // Создаем массив указателей
    // А дальше работа как с обычным массивом. 
    for (int i = 0; i < n; i++)
    {
        a[i] =  (double)(rand()) / RAND_MAX * 200 - 100;
        matrix << a[i] << " "; // Вывести элементы на консольку
    }
    delete [] a; // А потом массив
    matrix.close();
    return 0;
}