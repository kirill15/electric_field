/* Класс MatrixFEM - конечноэлементраная матрица
 */

#ifndef MATRIXFEM_H
#define MATRIXFEM_H

#include <set>
#include "grid2D.h"

// Тип - Матрица в разреженном строчном формате
/* Индексы с 0! */
struct Matrix
{
    unsigned int n;   // Размерность матрицы
    unsigned int *ig; // Указатели начала строк
    unsigned int *jg; // Номера столбоцов внедиагональных элементов матрицы
    double *ggl; 	  // Элементы нижнего треуголника
    double *ggu; 	  // Элементы верхнего треуголника
    double *di;		  // Диагональные элементы матрицы

    Matrix() : ig(nullptr), jg(nullptr), ggl(nullptr), ggu(nullptr), di(nullptr) {}

    ~Matrix()
    {/*
        delete[] ig;
        delete[] jg;
        delete[] ggl;
        //if (ggu)
          //  delete[] ggu;
        delete[] di;*/
    }
};


class MatrixFEM
{
private:
    Matrix a; // Матрица в разреженном строчно-столбцовом формате
    size_t size_jg; // Размер массива jg

public:
    MatrixFEM();


    // Выделяет память под элементы матрицы и обнуляет все элементы
    void createArrays();


    // Построение портрета матрицы СЛАУ
    /* Создает массивы ig и jg
     *
     * grid::Grid* - расчетная сетка
     * sizeOfElement::int - вид К.Э. (3 - треугольники, 4 - четырехугольники)
     */
    void generatePortrait(Grid2D *grid, unsigned sizeOfElement = 4);

    // Считывание портрета матрицы СЛАУ
    /* В файле ig первая запись - размерность матрицы
     */
    int readPortrait(string ig, string jg);

    // Сохранение портрета матрицы СЛАУ (ig.txt, jg.txt)
    void savePortrait();


    // Сохранение массивов ggl, ggu, di
    void saveElements();


    // Возвращает матрицу
    Matrix *matrix();
};

#endif // MATRIXFEM_H
