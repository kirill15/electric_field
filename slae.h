/* Класс SLAE будет содержать различные методы решения СЛАУ
 * ЛОС для несимметричной матрицы
 * МСГ для симметричой матрицы
*/

#ifndef SLAE_H
#define SLAE_H

#include <iostream>
#include "matrix.h"
#include <cmath>


class SLAE
{
private:
    SLAE();

    // Обновление вектора невязки
    static void refreshRk(Matrix &a, double *f, double *x, double *r);


public:
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

    // Норма вектора
    static double norm(const double *a, int size);


    // Неполная LU-факторизация матрицы
    static void factorizeLU(Matrix &a, Matrix &lu);

    // Неполная LL^T-факторизация матрицы
    static void factorizeLLT(Matrix &a, Matrix &l);


    // Прямой ход (Ly=b)
    static void forwardStroke(Matrix &lu, double *b, double *y);

    // Прямой ход для LLT (Ly=b)
    static void forwardStrokeLLT(Matrix &l, double *b, double *y);

    // Обратный ход (Ux=b)
    static void returnStroke(Matrix &lu, double *b, double *x);

    // Обратный ход для LLT (Ux=b)
    static void returnStrokeLLT(Matrix &l, double *b, double *x);


    // Итерационный шаг для матрицы c предобуславливанием неполным LU-разложением (ЛОС)
    static double iterationLU(Matrix &a, Matrix &lu, double *x, double *r, double *z, double *p);

    // ЛОС для матрицы c предобуславливанием неполным LU-разложением
    static int solveLOS_LU(Matrix &a, double *f, double *x, double eps, int maxIter = 1000);


    // Итерационный шаг для матрицы c предобуславливанием неполным LLT-разложением (МСГ)
    static void iterationMSGLLT(Matrix &a, Matrix &l, double *x, double *r, double *z, double *temp1, double *temp2);

    // МСГ для матрицы c предобуславливанием неполным LLT-разложением
    static int solveMSG_LLT(Matrix &a, double *f, double *x, double eps, int maxIter);
};

#endif // SLAE_H








