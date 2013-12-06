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

    double eps; // Точность решения СЛАУ (По-умолчанию 1e-15)

    unsigned sourceNode; // Номер узла точечного источника
    double sourceValue; // Значение плотности источника

    // Значение sigma
    double sigma(int nvk);


public:
    NormalField();

    ~NormalField()
    {
        delete grid;
        delete matrix;
        delete[] f;
        delete[] v;
    }


    // Локальная матрица жесткости
    /* p :: Coord[] - координаты вершин конечного элемента
     * nvk :: int - номер в каталоге
     * G :: double[][] - локальная матрица жесткости (возвращаемое значение)
     */
    void createLocalG(const Coord p[4], int nvk, double G[4][4]);

    // Глобальная Матрица
    void createGlobalMatrix();

    // Глобальный вектор правой части
    void createGlobalRightPart();

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

    // Получить значение eps
    double getEps() const;

    // Задать значение eps
    void setEps(double value);

    // Получить позицию и значение источника
    void getSource(unsigned &node, double &value) const;

    // Задать позицию и значение источника
    void setSource(unsigned node, double value);

    // Получить вектор весов решения
    double *getV(unsigned &size) const;
};

#endif // NORMALFIELD_H
