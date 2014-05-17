#ifndef ANOMALOUSFIELD_H
#define ANOMALOUSFIELD_H

#include "grid3D.h"
#include "matrixfem.h"
#include "normalfield.h"

using std::string;

class AnomalousField
{
private:
    Grid3D *grid; // Расчетная сетка

    MatrixFEM *matrix; // Конечноэлементная матрица

    double *f; // Правая часть

    Coord3D pSourcePlus; // Координаты источника (+)
    Coord3D pSourceMinus; // Координаты источника (-)

    NormalField *v0;

    double *v; // Вектор весов (решение СЛАУ)

    double eps; // Точность решения СЛАУ (По-умолчанию 1e-15)

    double **sigmas; // Массив значений сигм (индекс - номер подобласти)
    double **rotations; // Матрицы поворота для сигм

private:
    // Значение Ug |||||||||||||||||| Для теста |||||||||||||||||||
    /* rz::Coord - координата узла
     * nk::int - номер в каталоге первого краевого условия
     */
//    double GetUg(Coord3D p, int nvk);

    // V0(x, y, z)
    double getV0(Coord3D p);


    // Значение F |||||||||||||||||| Для теста |||||||||||||||||||
//    double getF(Coord3D p, int nvk);


    // Матрица жесткости одномерного элемента
    void localG_1D(double G[2][2], double h);

    // Матрица массы одномерного элемента
    void localM_1D(double M[2][2], double h);

    // Учет первого краевого на грани
    /* startIntex::size_t - стартовый узел
     * stepOutside::size_t - размер шага для перехода на следующий "уровень"
     * sizeInside::size_t - размер ребра
     * limit::size_t - максимальное значение номера узла
     * stepInside::size_t - шаг "внутренний"
     */
    void firstBoundaryConditionOnFace(size_t startIndex, size_t stepOutside, size_t sizeInside, size_t limit = 0, size_t stepInside = 1);


public:
    AnomalousField();


    // Локальная матрица жесткости
    /* p :: Coord3D[] - координаты вершин конечного элемента
     * nvk :: int - номер в каталоге
     * G :: double[][] - локальная матрица жесткости (возвращаемое значение)
     */
    void createLocalG(const Coord3D p[8], double *G[8][8]);

    // Локальный вектор правой части
    void createLocalF(const Coord3D p[8], int nvk, double *G[8][8], double F[8]);

    // Считывание значений сигмы (и ее матрицы поворота) из файла
    int readSigma(string fileWithSigma);


    // Глобальная СЛАУ
    void createGlobalSLAE();

    // Учет первого краевого условия
    void firstBoundaryCondition();

    // Решение СЛАУ
    /* method :: string - метод решения СЛАУ
     * maxIter :: size_t - максимальное число итераций
     * x0 :: double - начальное приближение (x0, x0, x0,...,x0)
     */
    void solve(string method = "MSG_LLT", size_t maxIter = 1000, double x0 = 0.0001);


    // Создание сетки из файлов, описывающих область и ее разбиение
    int createGrid(string fileWithArea, string fileWithGrid, size_t fragmentation = 0);

    // Создание портрета по сетке
    void createPortrait();

    // Получить вектор весов решения
    double *getV(unsigned &size) const;

    // Задать координаты ГЭЛа
    void setSource(const Coord3D &plus, const Coord3D &minus);

    // Задать основное поле
    void setNormalField(NormalField *v);

    // Решение в точке
    double getValue(Coord3D xyz);


    double getEps() const;
    void setEps(double value);
    Grid3D *getGrid() const;
    void setGrid(Grid3D *value);
    MatrixFEM *getMatrix() const;
    void setMatrix(MatrixFEM *value);
};

#endif // ANOMALOUSFIELD_H
