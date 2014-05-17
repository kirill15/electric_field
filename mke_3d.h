#ifndef ANOMALOUSFIELD_H
#define ANOMALOUSFIELD_H

#include "grid3D.h"
#include "matrixfem.h"

class Mke3D
{
private:
    Grid3D *grid; // Расчетная сетка

    MatrixFEM *matrix; // Конечноэлементная матрица

    double *f; // Правая часть

    double *v; // Вектор весов (решение СЛАУ)

    double eps; // Точность решения СЛАУ (По-умолчанию 1e-15)

    double *sigmas; // Массив значений сигм (индекс - номер подобласти)

private:
    // Значение Ug |||||||||||||||||| Для теста |||||||||||||||||||
    /* rz::Coord - координата узла
     * nk::int - номер в каталоге первого краевого условия
     */
    double GetUg(Coord3D p, int nvk);

    // Значение F |||||||||||||||||| Для теста |||||||||||||||||||
    double getF(Coord3D p, int nvk);


    // Матрица жесткости одномерного элемента
    void localG_1D(double G[2][2], double h);

    // Матрица массы одномерного элемента
    void localM_1D(double M[2][2], double h);

    // Учет первого краевого на гране
    /* startIntex::size_t - стартовый узел
     * stepOutside::size_t - размер шага для перехода на следующий "уровень"
     * sizeInside::size_t - размер ребра
     * limit::size_t - максимальное значение номера узла
     * stepInside::size_t - шаг "внутренний"
     */
    void firstBoundaryConditionOnFace(size_t startIndex, size_t stepOutside, size_t sizeInside, size_t limit = 0, size_t stepInside = 1);


public:
    Mke3D();


    // Локальная матрица жесткости
    /* p :: Coord3D[] - координаты вершин конечного элемента
     * nvk :: int - номер в каталоге
     * G :: double[][] - локальная матрица жесткости (возвращаемое значение)
     */
    void createLocalG(const Coord3D p[8], int nvk, double G[8][8]);

    // Локальный вектор правой части |||||||||||||||||| Для теста |||||||||||||||||||
    void createLocalF(const Coord3D p[8], int nvk, double F[8]);

    // Считывание значений сигмы на слоях из файла
    int readSigma(std::string fileWithSigma);


    // Глобальная СЛАУ
    void createGlobalSLAE();


    // Учет первого краевого условия
    void firstBoundaryCondition();

    // Решение СЛАУ
    /* Начальное приближение задается внутри (НЕ ОЧЕНЬ) */
    void solve();



    // Создание сетки из файлов, описывающих область и ее разбиение
    int createGrid(std::string fileWithArea, std::string fileWithGrid, size_t fragmentation = 0);

    // Создание портрета по сетке
    void createPortrait();

    // Решение в точке
    double getValue(Coord3D xyz);

    // Получить вектор весов решения
    double *getV(unsigned &size) const;

    double getEps() const;
    void setEps(double value);
    Grid3D *getGrid() const;
    void setGrid(Grid3D *value);
    MatrixFEM *getMatrix() const;
    void setMatrix(MatrixFEM *value);
};

#endif // ANOMALOUSFIELD_H
