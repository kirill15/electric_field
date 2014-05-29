/* Класс Grid2D создает 2D сетку (файлы nvtr, nvkat, rz) на основе расчетной области и разбиения ее на подобласти.
 * Позволяет сохранять файлы в бинарном формате для Тельмы
 */

#ifndef GREED_H
#define GREED_H

#include <fstream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <iostream>

// Тип - rz-координата точки
struct Coord
{
    Coord(){}
    Coord(double r, double z)
    {
        this->r = r;
        this->z = z;
    }
    double r;
    double z;
};

class Grid2D
{
private:
    int nRW; // Размерность массива rw
    int nZW; // Размерность массива zw
    double *rw; // Массив координат r всех границ подобластей
    double *zw; // Массив координат z всех границ подобластей
    unsigned l; // Число подобластей
    int **areas; // Области (пятерки целых чисел)

    std::vector < std::pair <int, double> > partitionsR; // Пары <Количество подынтервалов по оси r> <Коэффициент разрядки>
    std::vector < std::pair <int, double> > partitionsZ; // Пары <Количество подынтервалов по оси z> <Коэффициент разрядки>

    uint countPoints; // Количество узлов
    uint countFE; // Количество КЭ

    unsigned sizeR; // Количество узлов по оси r
    unsigned sizeZ; // Количество узлов по оси z

public:
    // Массивы, описывающие КЭ:
    Coord *rz; // Координаты узлов
    uint **nvtr; // Глобальные номера узлов каждого КЭ
    uint *nvkat; // Номер области для каждого КЭ



private:
    // Создание массивов nvtr, rz, nvkat
    void createArrays(double *r, double *z, long sizeR, long sizeZ);

public:
    Grid2D();
    virtual ~Grid2D();



    // Считывание расчетной области
    /* Формат:
     *  <Длина массива координат r всех границ подобластей (rw)>
     *  <Перечень элементов этого массива rw>
     *  <Длина массива координат z всех границ подобластей (zw)>
     *  <Перечень элементов этого массива zw>
     *  <L - количество подобластей>
     *  <Номер формул>
     *   <Номер элемента в массиве rw (Левая координата r)>
     *   <Номер элемента в массиве rw (Правая координата r)>
     *   <Номер элемента в массиве zw (Нижняя координата z)>
     *   <Номер элемента в массиве zw (Верхняя координата z)>
     *  ...
     *
     * Возвращает 0 в случае успешного выполнения
     */
    int readArea(std::string fileWithArea);



    // Считывание разбиения расчетной области (сетки)
    /* Формат:
     *  <Количество подынтервалов> <Коэффициент разрядки> ...
     *  ...
     *
     * Возвращает 0 в случае успешного выполнения
     */
    int readPartitions(std::string fileWithGrid);



    // Создание массивов nvtr, rz, nvkat
    void createGrid();



    // Печать сетки в файл
    /* Создает текстовые файлы:
     *  nvtr.txt
     *  nvkat.txt
     *  rz.txt
     */
    virtual void saveGrid();



    // Перевод файлов сетки в бинарный формат
    /* Конвертирует файлы следующим образом:
    *   nvtr.txt -> nvtr.dat
    *   nvkat.txt -> nvkat2d.dat
    *   rz.txt -> rz.dat
    *  Создает файл inf2tr.dat
    *
    * pathToProgram :: string - Путь к программе конвертации
    * Возвращает 0 в случае успешного выполнения
    */
    static int txtToDat(std::string pathToProgram);


    // Возвращает номер подобласти по координате (проверка только по z)
    size_t getArea(double z);


    // Возвращает количество узлов
    virtual unsigned getCountPoints() const;

    // Возвращает количество К.Э.
    virtual unsigned getCountFE() const;

    // Возвращает количество узлов по оси r
    unsigned getSizeR() const;

    // Возвращает количество узлов по оси z
    unsigned getSizeZ() const;

    // Возвращает nvtr
    virtual unsigned **getNvtr() const;
};


#endif // GREED_H
