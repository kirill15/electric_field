#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "grid2D.h"

using namespace std;



// Тип - xyz-координата точки
struct Coord3D
{
    Coord3D(){}
    Coord3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    double x;
    double y;
    double z;
};


class Grid3D : public Grid2D
{
private:
    int nXW; // Размерность массива xw
    int nYW; // Размерность массива yw
    int nZW; // Размерность массива zw
    double *xw; // Массив координат x всех границ подобластей
    double *yw; // Массив координат y всех границ подобластей
    double *zw; // Массив координат z всех границ подобластей
    unsigned l; // Число подобластей
    int **areas; // Области (семерки целых чисел)

    vector < pair <int, double> > partitionsX; // Пары <Количество подынтервалов по оси x> <Коэффициент разрядки>
    vector < pair <int, double> > partitionsY; // Пары <Количество подынтервалов по оси y> <Коэффициент разрядки>
    vector < pair <int, double> > partitionsZ; // Пары <Количество подынтервалов по оси z> <Коэффициент разрядки>

    uint countPoints; // Количество узлов
    uint countFE; // Количество КЭ

    unsigned sizeX; // Количество узлов по оси x
    unsigned sizeY; // Количество узлов по оси y
    unsigned sizeZ; // Количество узлов по оси z

    int compare (const void * a, const void * b)
    {
      return ( *(double*)a - *(double*)b );
    }

public:
    // Массивы, описывающие КЭ:
    Coord3D *xyz; // Координаты узлов
    uint **nvtr; // Глобальные номера узлов каждого КЭ
    uint *nvkat; // Номер области для каждого КЭ


private:
    // Создание массивов nvtr, xyz, nvkat
    void createArrays(double *x, double *y, double *z);

    // Получить номер подобласти, которой принадлежит К.Э.
    /* x: левая граница, правая граница
     * y: ближняя граница, дальняя граница
     * z: нижняя граница, верхняя граница
     */
    unsigned getAreaNumber(int xLeft, int xRight, int yNear, int yFar, int zLow, int zTop);
public:
    Grid3D();
    virtual ~Grid3D();



    // Считывание расчетной области
    /* Формат:
     *  <Длина массива координат x всех границ подобластей (xw)>
     *  <Перечень элементов этого массива xw>
     *  <Длина массива координат y всех границ подобластей (yw)>
     *  <Перечень элементов этого массива yw>
     *  <Длина массива координат z всех границ подобластей (zw)>
     *  <Перечень элементов этого массива zw>
     *  <L - количество подобластей>
     *  <Номер формул>
     *   <Номер элемента в массиве xw (Левая координата x)>
     *   <Номер элемента в массиве xw (Правая координата x)>
     *   <Номер элемента в массиве yw (Ближняя координата y)>
     *   <Номер элемента в массиве yw (Дальняя координата y)>
     *   <Номер элемента в массиве zw (Нижняя координата z)>
     *   <Номер элемента в массиве zw (Верхняя координата z)>
     *  ...
     *
     * Возвращает 0 в случае успешного выполнения
     */
    int readArea(string fileWithArea);



    // Считывание разбиения расчетной области (сетки)
    /* Формат:
     *  <Количество подынтервалов> <Коэффициент разрядки> ...
     *  ...
     *
     * Возвращает 0 в случае успешного выполнения
     */
    int readPartitions(string fileWithGrid, size_t fragmentation = 0);



    // Создание массивов nvtr, xyz, nvkat
    void createGrid();


    // Сохранение сетки в файлы
    void saveGrid();


    // Возвращает количество узлов
    unsigned getCountPoints() const;


    // Возвращает количество К.Э.
    unsigned getCountFE() const;


    unsigned **getNvtr() const;







    unsigned getSizeX() const;
    unsigned getSizeY() const;
    unsigned getSizeZ() const;
};


#endif // GRID3D_H
