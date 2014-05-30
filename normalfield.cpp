#include "normalfield.h"
#include "slae.h"

using namespace std;


void NormalField::gaussianElimination(size_t index, size_t *relatedNodes, size_t countRelatedNodes)
{
    static Matrix a = matrix->matrix();

    // Обнуляем столбец выше элемента и строку левее него
    for (size_t i = a.ig[index]; i < a.ig[index + 1]; i++)
        a.ggl[i] = 0.0;

    // Обнуляем правее и ниже
    for (size_t node = 0; node < countRelatedNodes; node++)
    {
        size_t indexInIg = a.ig[relatedNodes[node]];

        // Выбираем нужный столбец
        while (a.jg[indexInIg] != index)
            indexInIg++;

        a.ggl[indexInIg] = 0.0;
    }
}

NormalField::NormalField() : grid(nullptr), matrix(nullptr), f(nullptr), v(nullptr), eps(1e-15)
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
        cerr << "Error with opening sigma.txt" << endl;
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

    // Выделяем память для элементов матрицы (Там же они и обнуляются)
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
    Matrix a = matrix->matrix();
    unsigned sizeR = grid->getSizeR();

    size_t nodesRelated[3]; // Узлы, большие по номеру, чем "краевой" узел

    /// Задаем однородные краевые условия снизу
    nodesRelated[0] = sizeR;
    nodesRelated[1] = sizeR + 1;
    a.di[0] = 1.0;
    gaussianElimination(0, nodesRelated, 2);

    size_t endIndex = sizeR - 1;
    for (size_t index = 1; index < endIndex; index++)
    {
        a.di[index] = 1.0;

        nodesRelated[0] = index + sizeR - 1;
        nodesRelated[1] = nodesRelated[0] + 1;
        nodesRelated[2] = nodesRelated[1] + 1;

        gaussianElimination(index, nodesRelated, 3);
    }

    /// Задаем однородные краевые условия справа
    endIndex = a.n - 1;
    for (size_t index = sizeR - 1; index < endIndex; index += sizeR)
    {
        a.di[index] = 1.0;

        nodesRelated[0] = index + sizeR - 1;

        gaussianElimination(index, nodesRelated, 1);
    }
    a.di[endIndex] = 1.0;
    gaussianElimination(endIndex, nullptr, 0);




/*
    size_t startIndex = 0;
    size_t endIndex = sizeR;
    size_t step = 1;

    for (int i = 0; i < 2; i++)
    {
        for (size_t index = startIndex; index < endIndex; index += step)
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
        startIndex = 2 * sizeR - 1;
        endIndex = a.n;
        step = sizeR;
    }
*/




/*
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
*/

//    matrix->saveElements();


}


void NormalField::solve(string method, size_t maxIter)
{
    // Выделяем память под результат
    if (v)
        delete[] v;
    v = new double[matrix->matrix().n];

    // Начальное приближение
    for (size_t i = 0; i < matrix->matrix().n; i++)
        v[i] = 0.0001;

    // Решаем СЛАУ
    if (method == "MSG_LLT")
        SLAE::solveMSG_LLT(matrix->matrix(), f, v, eps, maxIter);
    else if (method == "LOS_LU")
        SLAE::solveLOS_LU(matrix->matrix(), f, v, eps, maxIter);

//    matrix->saveElements();
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

double NormalField::getSigma(double z) const
{
    size_t num = grid->getArea(z);
    return sigmas[num];
}


double NormalField::getValue(Coord rz) const
{
    Coord *coords = grid->rz;
    size_t sizeX = grid->getSizeR();
    size_t sizeY = grid->getSizeZ();

    if (rz.r < coords[0].r || rz.r > coords[sizeX - 1].r || rz.z < coords[0].z || rz.z > coords[sizeX * sizeY - 1].z)
        throw "Point does not lie in the field of";

    // Получаем номера координатных линий
    size_t p, s, i; // p - индекс координаты по Х, s - по Y, i - номер узла (на левой границе)
    size_t first = 1, last = sizeX; // last - номер элемента в массиве, СЛЕДУЮЩЕГО ЗА последним
    while (first < last)
    {
        size_t mid = first + (last - first) / 2;
        if (rz.r <= coords[mid].r)
            last = mid;
        else
            first = mid + 1;
    }
    p = last;

    first = 1, last = sizeY;
    while (first < last)
    {
        size_t mid = first + (last - first) / 2;
        if (rz.z <= coords[mid * sizeX].z)
            last = mid;
        else
            first = mid + 1;
    }
    i = last * sizeX;
    s = last;

    // Получаем глобальные номера базисных функций
    size_t k[4];
    k[1] = (s - 1) * sizeX + p;
    k[0] = k[1] - 1;
    k[3] = s * sizeX + p;
    k[2] = k[3] - 1;

    // Считаем значения базисных функций на КЭ
    double h = coords[p].r - coords[p - 1].r;
    double R1 = (coords[p].r - rz.r) / h;
    double R2 = (rz.r - coords[p - 1].r) / h;
    h = coords[i].z - coords[i - 1].z;
    double Z1 = (coords[i].z - rz.z) / h;
    double Z2 = (rz.z - coords[i - 1].z) / h;

    return v[k[0]] * R1*Z1  +  v[k[1]] * R2*Z1  +  v[k[2]] * R1*Z2  +  v[k[3]] * R2*Z2;
}


Grid2D *NormalField::getGrid() const
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
























/* //Поиск элемента перебором
double NormalField::getValue(Coord rz)
{
    Coord *coords = grid->rz;
    size_t sizeX = grid->getSizeR();
    size_t sizeY = grid->getSizeZ();

    if (rz.r < coords[0].r || rz.r > coords[sizeX - 1].r || rz.z < coords[0].z || rz.z > coords[sizeX * sizeY - 1].z)
        throw "Point does not lie in the field of";

    // Получаем номера координатных линий
    size_t p, s, i; // p - индекс координаты по Х, s - по Y, i - номер узла (на левой границе)
    for (p = 1; p < sizeX && rz.r > coords[p].r; p++);
    for (i = 0, s = 0; s < sizeY && rz.z > coords[i].z; i += sizeX, s++);

    ofstream file("ttLine", ios_base::app);
    file << p << " " << s << " " << i << endl;
    file.close();

    // Получаем глобальные номера базисных функций
    size_t k[4];
    k[1] = (s - 1) * sizeX + p;
    k[0] = k[1] - 1;
    k[3] = s * sizeX + p;
    k[2] = k[3] - 1;

    // Считаем значения базисных функций на КЭ
    double h = coords[p].r - coords[p - 1].r;
    double R1 = (coords[p].r - rz.r) / h;
    double R2 = (rz.r - coords[p - 1].r) / h;
    h = coords[i].z - coords[i - 1].z;
    double Z1 = (coords[i].z - rz.z) / h;
    double Z2 = (rz.z - coords[i - 1].z) / h;

    return v[k[0]] * R1*Z1  +  v[k[1]] * R2*Z1  +  v[k[2]] * R1*Z2  +  v[k[3]] * R2*Z2;
}
*/
