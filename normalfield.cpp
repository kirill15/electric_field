#include "normalfield.h"
#include "slae.h"


NormalField::NormalField() : grid(nullptr), matrix(nullptr), f(nullptr), v(nullptr), eps(1e-30)
{
}


void NormalField::getSource(unsigned &node, double &value) const
{
    node = sourceNode;
    value = sourceValue;
}

void NormalField::setSource(double value, unsigned node)
{
    // Положение точечного источника
    sourceNode = node ? node - 1 : matrix->matrix().n - grid->getSizeR();

    cout << "sourceNode: " << sourceNode + 1 << endl;

    // Его значение
    sourceValue = value;
}

double *NormalField::getV(unsigned &size) const
{
    size = matrix->matrix().n;
    return v;
}

int NormalField::readSigma(string fileWithSigma)
{
    ifstream file(fileWithSigma.c_str());
    if (!file.is_open())
    {
        cerr << "Error with opening sigma.txt";
        return -1;
    }

    int n, num;
    double val;
    file >> n;

    sigmas = new double[n];

    for (int i = 0; i < n; i++)
    {
        file >> num >> val;
        sigmas[num - 1] = val;
    }

    file.close();

    return 0;
}


unsigned int NormalField::getCountFE() const
{
    return grid->getCountFE();
}



void NormalField::createLocalG(const Coord p[], int nvk, double G[4][4])
{
    // Левое значение r
    double rp = p[0].r;

    // Шаги по осям z и r
    double hR = p[1].r - p[0].r;
    double hZ = p[2].z - p[0].z;

    double s = sigmas[nvk];

    double hrrp_hz = hR * rp / hZ;
    double hzrp_hr = hZ * rp / hR;
    double hrSq_hz = hR * hR / hZ;

    G[0][0] = s  *  (hrrp_hz / 3.0   +    hZ / 6.0    +   hzrp_hr / 3.0   +   hrSq_hz / 12.0);
    G[0][1] = s  *  (hrrp_hz / 6.0   -    hZ / 6.0    -   hzrp_hr / 3.0   +   hrSq_hz / 12.0);
    G[0][2] = s  *  (-hrrp_hz / 3.0  +    hZ / 12.0   +   hzrp_hr / 6.0   -   hrSq_hz / 12.0);
    G[0][3] = s  *  (-hrrp_hz / 6.0  -    hZ / 12.0   -   hzrp_hr / 6.0   -   hrSq_hz / 12.0);

    G[1][1] = s  *  (hrrp_hz / 3.0   +    hZ / 6.0    +   hzrp_hr / 3.0   +   hrSq_hz / 4.0);
    G[1][2] = s  *  (-hrrp_hz / 6.0  -    hZ / 12.0   -   hzrp_hr / 6.0   -   hrSq_hz / 12.0);
    G[1][3] = s  *  (-hrrp_hz / 3.0  +    hZ / 12.0   +   hzrp_hr / 6.0   -   hrSq_hz / 4.0);

    G[2][2] = s  *  (hrrp_hz / 3.0   +    hZ / 6.0    +   hzrp_hr / 3.0   +   hrSq_hz / 12.0);
    G[2][3] = s  *  (hrrp_hz / 6.0   -    hZ / 6.0    -   hzrp_hr / 3.0   +   hrSq_hz / 12.0);

    G[3][3] = s  *  (hrrp_hz / 3.0   +    hZ / 6.0    +   hzrp_hr / 3.0   +   hrSq_hz / 4.0);
}

void NormalField::createGlobalMatrix()
{
    double G[4][4]; // Локальная матрица жесткости

    // Сетка
    Coord *rz = grid->rz;
    unsigned **nvtr = grid->nvtr;
    unsigned *nvkat = grid->nvkat;

    // Выделяем память для элементов матрицы
    matrix->createArrays();

    // Матрица
    Matrix a = matrix->matrix();

    // Количество узлов К.Э.
    unsigned countFE = grid->getCountFE();

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
                unsigned left = a.ig[nvtr[k][j]]; // Левая граница поиска

                while (a.jg[left] != column)
                    left++;

                a.ggl[left] += G[i][j];
            }
        }
    }
}

void NormalField::createGlobalRightPart()
{
    size_t n = matrix->matrix().n;
    if (f)
        delete[] f;
    f = new double[n];
    for (size_t i = 0; i < n; i++)
        f[i] = 0;

    // Учитываем точечный источник
    f[sourceNode] = sourceValue;
}

void NormalField::firstBoundaryCondition()
{
    Matrix a =  matrix->matrix();

    /// Задаем однородные краевые условия снизу
    unsigned sizeR = grid->getSizeR();

    for (size_t index = 0; index < sizeR; index++)
    {
        a.di[index] = 1.0;
        f[index] = 0.0;

        // Обнуляем столбец выше элемента и строку левее него
        for (size_t i = a.ig[index]; i < a.ig[index + 1]; i++)
            a.ggl[i] = 0.0;

        // Обнуляем правее и ниже
        for (size_t i = index + 1; i < a.n; i++)
        {
            size_t left = a.ig[i]; // Начальное значение левой границы поиска
            size_t right = a.ig[i + 1] - 1; // Начальное значение правой границы поиска

            // Ищем нужный столбец
            while (a.jg[left] < index && left <= right)
                left++;

            // Если в этом столбце есть элемент, то присваиваем ему значение
            if (a.jg[left] == index)
                a.ggl[left] = 0.0;
        }
    }

    /// Задаем однородные краевые условия справа
    for (size_t index = 2 * sizeR - 1; index < a.n; index += sizeR)
    {
        a.di[index] = 1.0;
        f[index] = 0.0;

        // Обнуляем столбец выше элемента и строку левее него
        for (size_t i = a.ig[index]; i < a.ig[index + 1]; i++)
            a.ggl[i] = 0.0;

        // Обнуляем правее и ниже
        for (size_t i = index + 1; i < a.n; i++)
        {
            size_t left = a.ig[i]; // Начальное значение левой границы поиска
            size_t right = a.ig[i + 1] - 1; // Начальное значение правой границы поиска

            // Ищем нужный столбец
            while (a.jg[left] < index && left <= right)
                left++;

            // Если в этом столбце есть элемент, то присваиваем ему значение
            if (a.jg[left] == index)
                a.ggl[left] = 0.0;
        }
    }
}


void NormalField::solve()
{
    // Выделяем память под результат
    if (v)
        delete[] v;
    v = new double[matrix->matrix().n];

    // Начальное приближение
    for (size_t i = 0; i < matrix->matrix().n; i++)
        v[i] = 0.0001;

    // Решаем СЛАУ
    //SLAE::solveLOS_LU(matrix->matrix(), f, v, eps, 15000);
    SLAE::solveMSG_LLT(matrix->matrix(), f, v, eps, 1000);

    matrix->saveElements();
}

void NormalField::saveSolve()
{
    ofstream f;
    f.open("v.txt");
    size_t n = matrix->matrix().n;
    for (uint i = 0; i < n; i++)
        f << v[i] << endl;
    f.close();
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

double NormalField::getEps() const
{
    return eps;
}

void NormalField::setEps(double value)
{
    eps = value;
}



