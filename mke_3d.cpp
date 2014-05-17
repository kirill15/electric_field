#include "mke_3d.h"

using namespace std;


Mke3D::Mke3D() : grid(nullptr), matrix(nullptr), f(nullptr), v(nullptr), eps(1e-15)
{
}


double Mke3D::getEps() const
{
    return eps;
}

void Mke3D::setEps(double value)
{
    eps = value;
}

Grid3D *Mke3D::getGrid() const
{
    return grid;
}

void Mke3D::setGrid(Grid3D *value)
{
    grid = value;
}

MatrixFEM *Mke3D::getMatrix() const
{
    return matrix;
}

void Mke3D::setMatrix(MatrixFEM *value)
{
    matrix = value;
}
double Mke3D::GetUg(Coord3D p, int nvk)
{
    return p.x*p.x*p.x*p.x * p.y*p.y*p.y*p.y * p.z*p.z*p.z*p.z;
    //return p.x*p.x*p.x * p.y*p.y*p.y * p.z*p.z*p.z;
    //return p.x * p.x  *  p.y *p.y  *  p.z * p.z;
    //return exp(p.x * p.y * p.z);
    //return p.x * p.y * p.z;
}

double Mke3D::getF(Coord3D p, int nvk)
{
    return -12.0*p.x*p.x*p.y*p.y*p.z*p.z * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    //return -6.0*p.x*p.y*p.z * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    //return -2.0 * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    //return -2.0 * exp(p.x * p.y * p.z) * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    //return 0.0;
}

void Mke3D::localG_1D(double G[2][2], double h)
{
    G[0][0] = G[1][1] = 1.0 / h;
    G[0][1] = G[1][0] = -G[0][0];
}

void Mke3D::localM_1D(double M[2][2], double h)
{
    M[0][0] = M[1][1] = h / 3.0;
    M[0][1] = M[1][0] = h / 6.0;
}

void Mke3D::firstBoundaryConditionOnFace(size_t startIndex, size_t stepOutside, size_t sizeInside, size_t limit, size_t stepInside)
{
    Matrix a = matrix->matrix();

    if (!limit)
        limit = a.n;

    Coord3D *xyz = grid->xyz;
    uint *nvkat = grid->nvkat;

    for (size_t i = startIndex; i < limit; i += stepOutside)
        for (size_t j = 0; j < sizeInside; j += stepInside)
        {
            size_t el = i + j;
            cout << el + 1 << endl;

            a.di[el] = 1.0;
            double Ug = GetUg(xyz[el], nvkat[el]);
            f[el] = Ug;

            // Обнуляем столбец выше элемента и строку левее него
            for(size_t i = a.ig[el]; i < a.ig[el+1]; i++)
            {
                f[a.jg[i]] -= a.ggl[i] * Ug;	// Отнимаем от bi Aik*Ug
                a.ggl[i] = 0.0;                 // Обнуляем Aik
            }
            // Обнуляем правее и ниже
            for(size_t i = el+1; i < a.n; i++)
            {
                size_t left = a.ig[i],		// Начальное значение левой границы поиска
                       right = a.ig[i+1] - 1;	// Начальное значение правой границы поиска
                while(a.jg[left] < el && left <= right)	// Бинарный поиск нужного столбца (если он есть)
                {
                    size_t mid = (left + right) / 2;
                    if(a.jg[mid] < el)
                        left = mid + 1;
                    else
                        right = mid;
                }
                if(a.jg[left] == el)	// Проверяем, а вдруг в этом столбце вообще нет элемента
                {
                    f[i] -= a.ggl[left] * Ug;
                    a.ggl[left] = 0.0;
                }
            }
        }
    cout << endl;
}

int Mke3D::createGrid(string fileWithArea, string fileWithGrid, size_t fragmentation)
{
    if (grid)
        delete grid;
    else
        grid = new Grid3D;

    // Считываем область
    if (grid->readArea(fileWithArea))
    {
        cerr << "Не удалось считать область" << endl;
        return -1;
    }

    // Считываем разбиения области
    if (grid->readPartitions(fileWithGrid, fragmentation))
    {
        cerr << "Не удалось считать разбиение области" << endl;
        return -1;
    }

    // Создаем сетку
    grid->createGrid();

    return 0;
}

void Mke3D::createPortrait()
{
    if (matrix)
        delete matrix;
    else
        matrix = new MatrixFEM;

    matrix->generatePortrait(grid, 8);
}

double Mke3D::getValue(Coord3D xyz)
{
    Coord3D *coords = grid->xyz;
    size_t sizeX = grid->getSizeX();
    size_t sizeY = grid->getSizeY();
    size_t sizeZ = grid->getSizeZ();

    if (xyz.x < coords[0].x || xyz.x > coords[sizeX - 1].x || xyz.y < coords[0].y || xyz.y > coords[sizeX * sizeY - 1].y || xyz.z < coords[0].z || xyz.z > coords[sizeX * sizeY * sizeZ - 1].z)
        throw "Point does not lie in the field of";

    size_t sizeXY = sizeX * sizeY;

    // Получаем номера координатных линий
    size_t p, s, l, i; /* p - индекс координаты по Х, s - по Y, v - по Z */
    for (p = 1; p < sizeX && xyz.x > coords[p].x; p++);
    for (s = 0, i = 0; i < sizeY && xyz.y > coords[s].y; s += sizeX, i++);
    for (l = 0, i = 0; i < sizeZ && xyz.z > coords[l].z; l += sizeXY, i++);

    // Получаем глобальные номера базисных функций
    size_t k[8];
    k[1] = (l - sizeXY) + (s - sizeX) + p;
    k[0] = k[1] - 1;
    k[3] = k[1] + sizeX;
    k[2] = k[3] - 1;

    k[4] = k[0] + sizeXY;
    k[5] = k[1] + sizeXY;
    k[6] = k[2] + sizeXY;
    k[7] = k[3] + sizeXY;

    // Считаем значения базисных функций на КЭ
    double h = coords[p].x - coords[p - 1].x;
    double X1 = (coords[p].x - xyz.x) / h;
    double X2 = (xyz.x - coords[p - 1].x) / h;

    h = coords[s].y - coords[s - 1].y;
    double Y1 = (coords[s].y - xyz.y) / h;
    double Y2 = (xyz.y - coords[s - 1].y) / h;

    h = coords[l].z - coords[l - 1].z;
    double Z1 = (coords[l].z - xyz.z) / h;
    double Z2 = (xyz.z - coords[l - 1].z) / h;

    return v[k[0]] * X1*Y1*Z1  +  v[k[1]] * X2*Y1*Z1  +  v[k[2]] * X1*Y2*Z1  +  v[k[3]] * X2*Y2*Z1  +
           v[k[4]] * X1*Y1*Z2  +  v[k[5]] * X2*Y1*Z2  +  v[k[6]] * X1*Y2*Z2  +  v[k[7]] * X2*Y2*Z2;
}

double *Mke3D::getV(unsigned &size) const
{
    size = matrix->matrix().n;
    return v;
}

void Mke3D::createLocalG(const Coord3D p[8], int nvk, double G[8][8])
{
    double hX = p[1].x - p[0].x;
    double hY = p[2].y - p[0].y;
    double hZ = p[4].z - p[0].z;

    double sigma = sigmas[nvk];

    // Одномерные матрицы жесткости и массы
    double Gx[2][2];
    double Gy[2][2];
    double Gz[2][2];
    double Mx[2][2];
    double My[2][2];
    double Mz[2][2];

    localG_1D(Gx, hX);
    localG_1D(Gy, hY);
    localG_1D(Gz, hZ);
    localM_1D(Mx, hX);
    localM_1D(My, hY);
    localM_1D(Mz, hZ);


    size_t mui, nui, vi, muj, nuj, vj;

    for (size_t i = 0; i < 8; i++)
    {
        mui = i % 2;
        nui = (i / 2) % 2;
        vi = (i / 4);
        for (size_t j = 0; j < 8; j++)
        {
            muj = j % 2;
            nuj = (j / 2) % 2;
            vj = (j / 4);

            G[i][j] = sigma * (Gx[mui][muj] * My[nui][nuj] * Mz[vi][vj] +
                               Mx[mui][muj] * Gy[nui][nuj] * Mz[vi][vj] +
                               Mx[mui][muj] * My[nui][nuj] * Gz[vi][vj]);
        }
    }
}

void Mke3D::createLocalF(const Coord3D p[], int nvk, double F[])
{
    double hX = p[1].x - p[0].x;
    double hY = p[2].y - p[0].y;
    double hZ = p[4].z - p[0].z;

    // Одномерные матрицы массы
    double Mx[2][2];
    double My[2][2];
    double Mz[2][2];

    localM_1D(Mx, hX);
    localM_1D(My, hY);
    localM_1D(Mz, hZ);

    // Трехмерная матрица массы
    double C[8][8];
    size_t mui, nui, vi, muj, nuj, vj;
    for (size_t i = 0; i < 8; i++)
    {
        mui = i % 2;
        nui = (i / 2) % 2;
        vi = (i / 4);
        for (size_t j = 0; j < 8; j++)
        {
            muj = j % 2;
            nuj = (j / 2) % 2;
            vj = (j / 4);

            C[i][j] = Mx[mui][muj] * My[nui][nuj] * Mz[vi][vj];
        }
    }

    // f
    double f[8];
    for (size_t i = 0; i < 8; i++)
        f[i] = getF(p[i], nvk);

    // F=M*f
    for (size_t i = 0; i < 8; i++)
    {
        double s = 0.0;
        for (size_t j = 0; j < 8; j++)
            s += C[i][j] * f[j];
        F[i] = s;
    }
}

int Mke3D::readSigma(string fileWithSigma)
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

void Mke3D::createGlobalSLAE()
{
    double G[8][8]; // Локальная матрица жесткости
    double F[8]; // Локальный вектор правой части

    // Сетка
    Coord3D *xyz = grid->xyz;
    unsigned **nvtr = grid->nvtr;
    unsigned *nvkat = grid->nvkat;

    // Выделяем память для элементов матрицы (Там же они и обнуляются)
    matrix->createArrays();

    // Матрица
    Matrix a = matrix->matrix();

    // Память для вектора правой части
    f = new double[a.n];
    for (size_t i = 0; i < a.n; i++)
        f[i] = 0.0;

    // Количество узлов К.Э.
    unsigned countFE = grid->getCountFE();

    Coord3D p[8]; // Узлы К.Э.

    // Цикл по конечным элементам
    for (size_t k = 0; k < countFE; k++)
    {
        // Узлы конечного элемента
        for (size_t i = 0; i < 8; i++)
            p[i] = xyz[nvtr[k][i]];

        // Создаем локальную матрицу жесткости
        createLocalG(p, nvkat[k], G);

        // Создаем локальный вектор правой части
        createLocalF(p, nvkat[k], F);

        // Заносим матрицу жесткости в глобальную СЛАУ
        for(int i = 0; i < 8; i++)	// Цикл по строкам локальной матрицы
        {
            a.di[nvtr[k][i]] += G[i][i];    // Диагональные элементы

            uint column = nvtr[k][i];		// Искомый номер столбца

            for(int j = i + 1; j < 8; j++)	// Цикл по столбцам (верхнего треугольника)
            {
                uint left = a.ig[nvtr[k][j]];		// Начальное значение левой границы поиска
                uint right = a.ig[nvtr[k][j]+1] - 1;	// Начальное значение правой границы поиска

                while(a.jg[left] != column)			// Бинарный поиск нужного столбца
                {
                    int mid = (left + right) / 2;
                    if(a.jg[mid] < column)
                        left = mid + 1;
                    else
                        right = mid;
                }
                a.ggl[left] += G[i][j];
            }
        }

        // Заносим локальный вектор правой части в глобальный
        for (size_t i = 0; i < 8; i++)
            f[nvtr[k][i]] += F[i];
    }

    cout << "f:\n";
    for (size_t i = 0; i < a.n; i++)
        cout << f[i] << endl;
}

void Mke3D::firstBoundaryCondition()
{
    Matrix a = matrix->matrix();

    size_t sizeX = grid->getSizeX();
    size_t sizeY = grid->getSizeY();

    size_t sizeXY = sizeX * sizeY;

    // Ближняя грань
    firstBoundaryConditionOnFace(0, sizeXY, sizeX);
    // Дальняя грань
    firstBoundaryConditionOnFace(sizeXY - sizeX, sizeXY, sizeX);
    // Левая грань
    firstBoundaryConditionOnFace(sizeX, sizeXY, sizeXY - 2 * sizeX, a.n, sizeX);
    // Правая грань
    firstBoundaryConditionOnFace(2 * sizeX - 1, sizeXY, sizeXY - 2 * sizeX, a.n, sizeX);
    // Нижняя грань
    firstBoundaryConditionOnFace(sizeX + 1, sizeX, sizeX - 2, sizeXY - sizeX - 1);
    // Верхняя грань
    firstBoundaryConditionOnFace(a.n - sizeXY + sizeX + 1, sizeX, sizeX - 2, a.n - sizeX);
}

void Mke3D::solve()
{
    // Выделяем память под результат
    if (v)
        delete[] v;
    v = new double[matrix->matrix().n];

    // Начальное приближение
    for (size_t i = 0; i < matrix->matrix().n; i++)
        v[i] = 0.0001;

    // Решаем СЛАУ
    //SLAE::solveLOS_LU(matrix->matrix(), f, v, eps, 1000);
    SLAE::solveMSG_LLT(matrix->matrix(), f, v, eps, 1000);
}

















