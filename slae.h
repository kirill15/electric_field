/* Класс SLAE будет содержать различные методы решения СЛАУ
 * Пока содержит ЛОС для несимметричной матрицы
*/

#ifndef SLAE_H
#define SLAE_H

#include <iostream>
#include "matrix.h"

using namespace std;


class SLAE
{
private:


public:
    SLAE();

    // Умножение матрицы на вектор
    static void multMatrixVector(Matrix &a, double *x, double *res);

    // Сложение векторов (res = a + b)
    static void addVectorVector(double *a, double *b, int size, double *res);

    // Вычитание векторов (res = a - b)
    static void subVectorVector(double *a, double *b, int size, double *res);

    // Умножение вектора на число
    static void multVectorScalar(const double *v, double c, int size, double *res);

    // Скалярное произведение (в евклидовом пространстве)
    static double scalarProduct(const double *a, const double *b, int size);


    // Неполная LU-факторизация матрицы
    static void factorizeLU(Matrix &a, Matrix &lu);


    // Прямой ход (Ly=b)
    static void forwardStroke(Matrix &lu, double *b, double *y);

    // Обратный ход (Ux=b)
    static void returnStroke(Matrix &lu, double *b, double *x);


    // Итерационный шаг для матрицы c предобуславливанием неполным LU-разложением (ЛОС)
    static double iterationLU(Matrix &a, Matrix &lu, double *x, double *r, double *z, double *p);

    // ЛОС для матрицы c предобуславливанием неполным LU-разложением
    static int solveLOS_LU(Matrix &a, double *f, double *x, double eps, int maxIter = 1000);
};

#endif // SLAE_H








