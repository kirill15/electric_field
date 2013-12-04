#include "matrixfem.h"

MatrixFEM::MatrixFEM()
{
}

void MatrixFEM::createArrays()
{
    a.di = new double[a.n];
    a.ggl = new double[size_jg];
    a.ggu = a.ggl;

    for (size_t i = 0; i < a.n; i++)
        a.di[i] = 0.0;

    for (size_t i = 0; i < size_jg; i++)
        a.ggl[i] = 0.0;
}


Matrix* MatrixFEM::matrix()
{
    return &a;
}


void MatrixFEM::generatePortrait(Grid2D *grid, unsigned sizeOfElement)
{
    // Размерность матрицы
    a.n = grid->getCountPoints();

    // Создание списка смежности для каждой вершины
    set<unsigned> *list = new set<unsigned>[a.n];

    // Список смежных узлов для каждого узла в порядке возрастания
    size_t countFE = grid->getCountFE();
    for (size_t i = 0; i < countFE; i++)
    {
        uint *tmp = grid->nvtr[i];

        // Для каждого узла текущего элемена...
        for (size_t j = 0; j < sizeOfElement; j++)
        {
            // ...добавляем в список смежные с ним узлы
            for (size_t k = 0; k < sizeOfElement; k++)
            {
                // Добавляем только те узлы, глобальных номер которых больше текущего
                if (tmp[k] < tmp[j])
                    list[tmp[j]].insert(tmp[k]);
            }
        }
    }

    // Создание массива указателей начала строк ig
    if (a.ig)
        delete[] a.ig;
    else
        a.ig = new unsigned int[a.n + 1];

    a.ig[0] = 0;
    for (size_t i = 0; i < a.n; i++)
        // Прибавляем количество смежных вершин
        a.ig[i + 1] = a.ig[i] + list[i].size();

    // Создание массива указателей начала столбцов jg
    if (a.jg)
        delete[] a.jg;
    else
        a.jg = new unsigned int[a.ig[a.n]];

    size_jg = 0;
    for (size_t i = 1; i < a.n; i++)
        // Проходим по смежным элементам и добавляем их в массив jg
        for (auto j : list[i])
            a.jg[size_jg++] = j;

    delete[] list;
}

int MatrixFEM::readPortrait(string ig, string jg)
{
    ifstream f(ig.c_str());
    if (!f.is_open())
        return -1;

    f >> a.n;

    if (a.ig)
        delete[] a.ig;
    a.ig = new unsigned[a.n + 1];
    for (size_t i = 0; i <= a.n; i++)
    {
        f >> a.ig[i];
        a.ig[i]--;
    }
    f.close();

    f.open(jg.c_str());
    if (!f.is_open())
        return -2;

    size_jg = a.ig[a.n];

    if (a.jg)
        delete[] a.jg;
    a.ig = new unsigned[size_jg];
    for (size_t i = 0; i < size_jg; i++)
    {
        f >> a.jg[i];
        a.jg[i]--;
    }
    f.close();

    return 0;
}


void MatrixFEM::savePortrait()
{
    ofstream f("ig.txt");

    f << a.n << endl << endl;
    for (size_t i = 0; i <= a.n; i++)
        f << a.ig[i] + 1 << endl;
    f.close();

    f.open("jg.txt");
    for (size_t i = 0; i < size_jg; i++)
        f << a.jg[i] + 1 << endl;
    f.close();
}

void MatrixFEM::saveElements()
{
    ofstream f;
    f.open("ggl.txt");

    for (size_t i = 0; i < size_jg; i++)
        f << a.ggl[i] << endl;
    f.close();

    f.open("ggu.txt");
    for (size_t i = 0; i < size_jg; i++)
        f << a.ggu[i] << endl;
    f.close();

    f.open("di.txt");
    for (size_t i = 0; i < a.n; i++)
        f << a.di[i] << endl;
    f.close();
}







