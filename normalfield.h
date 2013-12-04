#ifndef NORMALFIELD_H
#define NORMALFIELD_H

#include "grid2D.h"
#include "matrixfem.h"

class NormalField
{
private:
    unsigned int countPoints; // Количество узлов
    unsigned int countFE; // Количество КЭ

    Grid2D *grid; // Расчетная сетка

    MatrixFEM *matrix; // Конечноэлементная матрица

    double *f; // Правая часть

    double *v; // Вектор весов (решение СЛАУ)

    double eps;

    // Значение sigma
    double sigma(int nvk);

public:
    NormalField();

    ~NormalField()
    {
        delete grid;
        delete matrix;
    }


    // Локальная матрица жесткости
    /* p :: Coord[] - координаты вершин конечного элемента
     * nvk :: int - номер в каталоге
     * G :: double[][] - локальная матрица жесткости (возвращаемое значение)
     */
    void createLocalG(const Coord p[4], int nvk, double G[4][4]);

    // Глобальная СЛАУ
    void createGlobalSLAE();

    // Решение СЛАУ
    void solve();


    // Создание сетки из файлов, описывающих область и ее разбиение
    int createGrid(string fileWithArea, string fileWithGrid);

    // Создание сетки из файлов nvtr, nvkat, rz
    int createGrid(string nvtr, string nvkat, string rz);


    // Создание портрета по сетке
    void createPortrait();

    // Создание портрета из файлов ig, jg
    /******************************** РЕАЛИЗОВАТЬ! *********************************/
    int createPortrait(string ig, string jg);





    // Получить сетку
    Grid2D *getGrid();

    // Задать сетку
    void setGrid(Grid2D *value);


    // Получить матрицу
    MatrixFEM *getMatrix() const;

    // Задать матрицу
    void setMatrix(MatrixFEM *value);
};

#endif // NORMALFIELD_H
