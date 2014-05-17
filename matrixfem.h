/* Класс MatrixFEM - конечноэлементраная матрица
 */

#ifndef MATRIXFEM_H
#define MATRIXFEM_H

#include <set>
#include "grid2D.h"
#include "slae.h"

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
    /* В файле ig первая запись - размерность матрицы (****НЕ ХАРОШИЙ****)
     */
    int readPortrait(std::string ig, std::string jg);

    // Сохранение портрета матрицы СЛАУ (ig.txt, jg.txt)
    void savePortrait();


    // Сохранение массивов ggl, ggu, di
    void saveElements();


    // Возвращает матрицу
    Matrix &matrix();
};

#endif // MATRIXFEM_H
