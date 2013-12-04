#include "normalfield.h"
#include "slae.h"


double NormalField::sigma(int nvk)
{
    return 1.0;
}

NormalField::NormalField() : grid(nullptr), matrix(nullptr)
{
    eps = 1e-15;
}

void NormalField::createLocalG(const Coord p[], int nvk, double G[4][4])
{
    // Левое значение r
    double rp = p[0].r;

    // Шаги по осям z и r
    double hR = p[1].r - p[0].r;
    double hZ = p[2].z - p[0].z;

    double s = sigma(nvk);

    G[0][0] = s  *  ((hR * rp) / (3.0 * hZ)   +    hZ / 6.0    +   (hZ * rp) / (3.0 * hR)   +   (hR * hR) / (12.0 * hZ));
    G[0][1] = s  *  ((hR * rp) / (6.0 * hZ)   -    hZ / 6.0    -   (hZ * rp) / (3.0 * hR)   +   (hR * hR) / (12.0 * hZ));
    G[0][2] = s  *  (-(hR * rp) / (3.0 * hZ)   +   hZ / 12.0   +   (hZ * rp) / (6.0 * hR)   -   (hR * hR) / (12.0 * hZ));
    G[0][3] = s  *  (-(hR * rp) / (6.0 * hZ)   -   hZ / 12.0   -   (hZ * rp) / (6.0 * hR)   -   (hR * hR) / (12.0 * hZ));

    G[1][0] = G[0][1];
    G[1][1] = s  *  ((hR * rp) / (3.0 * hZ)   +   hZ / 6.0   +   (hZ * rp) / (3.0 * hR)   +   (hR * hR) / (4.0 * hZ));
    G[1][2] = s  *  (-(hR * rp) / (6.0 * hZ)   -   hZ / 12.0   -   (hZ * rp) / (6.0 * hR)   -   (hR * hR) / (12.0 * hZ));
    G[1][3] = s  *  (-(hR * rp) / (3.0 * hZ)   +   hZ / 12.0   +   (hZ * rp) / (6.0 * hR)   -   (hR * hR) / (4.0 * hZ));

    G[2][0] = G[0][2];
    G[2][1] = G[1][2];
    G[2][2] = s  *  ((hR * rp) / (3.0 * hZ)   +   hZ / 6.0   +   (hZ * rp) / (3.0 * hR)   +   (hR * hR) / (12.0 * hZ));
    G[2][3] = s  *  ((hR * rp) / (6.0 * hZ)   -   hZ / 6.0   -   (hZ * rp) / (3.0 * hR)   +   (hR * hR) / (12.0 * hZ));

    G[3][0] = G[0][3];
    G[3][1] = G[1][3];
    G[3][2] = G[2][3];
    G[3][3] = s  *  ((hR * rp) / (3.0 * hZ)   +   hZ / 6.0   +   (hZ * rp) / (3.0 * hR)   +   (hR * hR) / (4.0 * hZ));
}

void NormalField::createGlobalSLAE()
{
    double G[4][4]; // Локальная матрица жесткости

    // Сетка
    Coord *rz = grid->rz;
    unsigned **nvtr = grid->nvtr;
    unsigned *nvkat = grid->nvkat;

    // Выделяем память для элементов матрицы
    matrix->createArrays();

    // Матрица
    Matrix a = *(matrix->matrix());

    Coord p[4]; // Узлы К.Э.

    // Цикл по конечным элементам
    for (size_t k = 0; k < countFE; k++)
    {
        // Узлы конечного элемента
        for (size_t i = 0; i < 4; i++)
            p[i] = rz[nvtr[k][i]];

        // Создаем локальную матрицу жесткости
        createLocalG(p, nvkat[k], G);

        // Заносим матрицу жесткости в глобальную СЛАУ
        for (size_t i = 0; i < 4; i++) // Цикл по строкам локальной матрицы
        {
            // Диагональные элементы
            a.di[nvtr[k][i]] += G[i][i];

            // Искомый номер столбца
            unsigned column = nvtr[k][i];

            //**** Тупо, да и ваще строчно-столбцовый формат хранения матрицы в данном случае не нужен ****//
            for (size_t j = i + 1; j < 4; j++) // Цикл по столбцам локальной матрицы
            {
                unsigned positionInGgl = a.ig[nvtr[k][j]];

                while (a.jg[positionInGgl] != column)
                    positionInGgl++;

                a.ggl[positionInGgl] += G[i][j];
            }
        }
    }

    // Правая часть
    b = new double[a.n];
    for (size_t i = 0; i < a.n; i++)
        b[i] = 0.0;
}

void NormalField::solve()
{
    SLAE::solveLOS_LU(matrix->matrix(), f, v, eps);
}


int NormalField::createGrid(string fileWithArea, string fileWithGrid)
{
    if (grid)
        delete grid;
    else
        grid = new Grid2D;

    // Считываем область
    if (grid->readArea(fileWithArea))
    {
        cerr << "Не удалось считать область" << endl;
        return -1;
    }

    // Считываем разбиения области
    if (grid->readPartitions(fileWithGrid))
    {
        cerr << "Не удалось считать разбиение области" << endl;
        return -1;
    }

    // Создаем сетку
    grid->createGrid();

    // Сохраняем количество узлов и К.Э. (************ ЗАЧЕМ? ***********)
    countPoints = grid->getCountPoints();
    countFE = grid->getCountFE();

    return 0;
}


void NormalField::createPortrait()
{
    if (matrix)
        delete matrix;
    else
        matrix = new MatrixFEM;

    matrix->generatePortrait(grid, 4);
}




Grid2D *NormalField::getGrid()
{
    return grid;
}

void NormalField::setGrid(Grid2D *value)
{
    if (grid)
        delete grid;

    grid = value;
}

MatrixFEM *NormalField::getMatrix() const
{
    return matrix;
}

void NormalField::setMatrix(MatrixFEM *value)
{
    if (matrix)
        delete matrix;

    matrix = value;
}






