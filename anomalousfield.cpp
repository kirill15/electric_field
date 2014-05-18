#include "anomalousfield.h"

using namespace std;


AnomalousField::AnomalousField() : grid(nullptr), matrix(nullptr), f(nullptr), v0(nullptr), v(nullptr), eps(1e-15)
{
}


double AnomalousField::getEps() const
{
    return eps;
}

void AnomalousField::setEps(double value)
{
    eps = value;
}

Grid3D *AnomalousField::getGrid() const
{
    return grid;
}

void AnomalousField::setGrid(Grid3D *value)
{
    grid = value;
}

MatrixFEM *AnomalousField::getMatrix() const
{
    return matrix;
}

void AnomalousField::setMatrix(MatrixFEM *value)
{
    matrix = value;
}

/*double AnomalousField::GetUg(Coord3D p, int nvk)
{
    return 0.0;

    //return 1.0 / (2.0 * M_PI * sqrt(p.x * p.x + p.y *p.y + p.z * p.z) * 0.01);
    //return p.x * p.x * p.y *p.y * p.z * p.z;
    //return exp(p.x * p.y * p.z);
    //return p.x * p.y * p.z;
}*/

double AnomalousField::getV0(Coord3D p)
{
    double r = sqrt((p.x - pSourcePlus.x) * (p.x - pSourcePlus.x)  +  (p.y - pSourcePlus.y) * (p.y - pSourcePlus.y));
    Coord rz(r, p.z);
    double v0B = v0->getValue(rz);

    rz.r = sqrt((p.x - pSourceMinus.x) * (p.x - pSourceMinus.x)  +  (p.y - pSourceMinus.y) * (p.y - pSourceMinus.y));
    double v0A = v0->getValue(rz);

    return v0B - v0A;
}

/*
double AnomalousField::getF(Coord3D p, int nvk)
{
    return 0.0;
    //return -2.0 * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    //return -2.0 * exp(p.x * p.y * p.z) * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
}
*/

void AnomalousField::localG_1D(double G[2][2], double h)
{
    G[0][0] = G[1][1] = 1.0 / h;
    G[0][1] = G[1][0] = -G[0][0];
}

void AnomalousField::localM_1D(double M[2][2], double h)
{
    M[0][0] = M[1][1] = h / 3.0;
    M[0][1] = M[1][0] = h / 6.0;
}

void AnomalousField::firstBoundaryConditionOnFace(size_t startIndex, size_t stepOutside, size_t sizeInside, size_t limit, size_t stepInside)
{
    Matrix a = matrix->matrix();

    if (!limit)
        limit = a.n;

    for (size_t i = startIndex; i < limit; i += stepOutside)
        for (size_t j = 0; j < sizeInside; j += stepInside)
        {
            size_t el = i + j;

            a.di[el] = 1.0;
            f[el] = 0.0;

            // Обнуляем столбец выше элемента и строку левее него
            for(size_t i = a.ig[el]; i < a.ig[el+1]; i++)
                a.ggl[i] = 0.0;                 // Обнуляем Aik

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
                    a.ggl[left] = 0.0;
            }
        }





    // Без учета однородности:
    /*
    Coord3D *xyz = grid->xyz;
    uint *nvkat = grid->nvkat;
    for (size_t i = startIndex; i < limit; i += stepOutside)
        for (size_t j = 0; j < sizeInside; j += stepInside)
        {
            size_t el = i + j;
//            cout << el + 1 << endl;

//            if (xyz[el] == Coord3D(0.0, 0.0, 0.0)) continue;

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
        }*/
}

int AnomalousField::createGrid(string fileWithArea, string fileWithGrid, size_t fragmentation)
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

void AnomalousField::createPortrait()
{
    if (matrix)
        delete matrix;
    else
        matrix = new MatrixFEM;

    matrix->generatePortrait(grid, 8);
}

double *AnomalousField::getV(unsigned &size) const
{
    size = matrix->matrix().n;
    return v;
}

void AnomalousField::setSource(const Coord3D &plus, const Coord3D &minus)
{
    pSourcePlus = plus;
    pSourceMinus = minus;
}

void AnomalousField::setNormalField(NormalField *v)
{
    v0 = v;
}

double AnomalousField::getValue(Coord3D xyz)
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
void AnomalousField::createLocalG(const Coord3D p[8], double *G[8][8])
{

    double hX = p[1].x - p[0].x;
    double hY = p[2].y - p[0].y;
    double hZ = p[4].z - p[0].z;

    double hyhz_hx = hY * hZ / (36.0 * hX);
    double hxhz_hy = hX * hZ / (36.0 * hY);
    double hxhy_hz = hX * hY / (36.0 * hZ);
    double hx_24 = hX / 24.0;
    double hy_24 = hY / 24.0;
    double hz_24 = hZ / 24.0;

    double K1[8][9] = {{ 4, 2, 2,   2, 4, 2,   2, 2, 4},
                       {-4, 2, 2,  -2, 2, 1,  -2, 1, 2},
                       { 2,-2, 1,   2,-4, 2,   1,-2, 2},
                       {-2,-2, 1,  -2,-2, 1,  -1,-1, 1},
                       { 2, 1,-2,   1, 2,-2,   2, 2,-4},
                       {-2, 1,-2,  -1, 1,-1,  -2, 1,-2},
                       { 1,-1,-1,   1,-2,-2,   1,-2,-2},
                       {-1,-1,-1,  -1,-1,-1,  -1,-1,-1}};

    double K2[7][9] = {{ 4,-2,-2,  -2, 4, 2,  -2, 2, 4},
                       {-2, 2,-1,   2,-2, 1,   1,-1, 1},
                       { 2, 2,-1,  -2,-4, 2,  -1,-2, 2},
                       {-2,-1, 2,   1, 1,-1,   2, 1,-2},
                       { 2,-1, 2,  -1, 2,-2,  -2, 2,-4},
                       {-1, 1, 1,   1,-1,-1,   1,-1,-1},
                       { 1, 1, 1,  -1,-2,-2,  -1,-2,-2}};

    double K3[6][9] = {{ 4,-2, 2,  -2, 4,-2,   2,-2, 4},
                       {-4,-2, 2,   2, 2,-1,  -2,-1, 2},
                       { 1, 1,-1,  -1,-2, 2,   1, 2,-2},
                       {-1, 1,-1,   1,-1, 1,  -1, 1,-1},
                       { 2,-1,-2,  -1, 2, 2,   2,-2,-4},
                       {-2,-1,-2,   1, 1, 1,  -2,-1,-2}};

    double K4[5][9] = {{ 4, 2,-2,   2, 4,-2,  -2,-2, 4},
                       {-1,-1, 1,  -1,-1, 1,   1, 1,-1},
                       { 1,-1, 1,   1,-2, 2,  -1, 2,-2},
                       {-2, 1, 2,  -1, 1, 1,   2,-1,-2},
                       { 2, 1, 2,   1, 2, 2,  -2,-2,-4}};

    double K5[4][9] = {{ 4, 2,-2,   2, 4,-2,  -2,-2, 4},
                       {-4, 2,-2,  -2, 2,-1,   2,-1, 2},
                       { 2,-2,-1,   2,-4,-2,  -1, 2, 2},
                       {-2,-2,-1,  -2,-2,-1,   1, 1, 1}};

    double K6[3][9] = {{ 4,-2, 2,  -2, 4,-2,   2,-2, 4},
                       {-2, 2, 1,   2,-2,-1,  -1, 1, 1},
                       { 2, 2, 1,  -2,-4,-2,   1, 2, 2}};

    double K7[2][9] = {{ 4,-2,-2,  -2, 4, 2,  -2, 2, 4},
                       {-4,-2,-2,   2, 2, 1,   2, 1, 2}};

    double K8[9]    =  { 4, 2, 2,   2, 4, 2,   2, 2, 4};

    double P[9];
    P[0] = hyhz_hx;
    P[1] = hz_24;
    P[2] = hy_24;
    P[3] = hz_24;
    P[4] = hxhz_hy;
    P[5] = hx_24;
    P[6] = hy_24;
    P[7] = hx_24;
    P[8] = hxhy_hz;

    for (size_t i = 0; i < 8; i++)
        for (size_t j = 0; j < 9; j++)
            G[0][i][j] = K1[i][j] * P[j];

    for (size_t i = 0; i < 7; i++)
        for (size_t j = 0; j < 9; j++)
            G[1][i+1][j] = K2[i][j] * P[j];

    for (size_t i = 0; i < 6; i++)
        for (size_t j = 0; j < 9; j++)
            G[2][i+2][j] = K3[i][j] * P[j];

    for (size_t i = 0; i < 5; i++)
        for (size_t j = 0; j < 9; j++)
            G[3][i+3][j] = K4[i][j] * P[j];

    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 9; j++)
            G[4][i+4][j] = K5[i][j] * P[j];

    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 9; j++)
            G[5][i+5][j] = K6[i][j] * P[j];

    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 9; j++)
            G[6][i+6][j] = K7[i][j] * P[j];

    for (size_t j = 0; j < 9; j++)
        G[7][7][j] = K8[j] * P[j];

}

void AnomalousField::createLocalF(const Coord3D p[8], int nvk, double *G[8][8], double F[])
{
    double *sigma = sigmas[nvk];
    double sigma0 = v0->getSigma((p[4].z + p[0].z) / 2.0);
    double dSigma[9];

    for (size_t i = 0; i < 9; i++)
        dSigma[i] = -sigma[i];
    for (size_t i = 0; i < 9; i++)
        dSigma[i] += rotations[nvk][i] * sigma0;

    double Gd[8][8];
    for (size_t i = 0; i < 8; i++)
        for (size_t j = i; j < 8; j++)
        {
            Gd[i][j] = 0.0;
            for (size_t l = 0; l < 9; l++)
                Gd[i][j] += G[i][j][l] * dSigma[l];
            Gd[j][i] = Gd[i][j];
        }

    // f
    double q0[8];
    for (size_t i = 0; i < 8; i++)
        q0[i] = getV0(p[i]);

    // F=G*q0
    for (size_t i = 0; i < 8; i++)
    {
        double s = 0.0;
        for (size_t j = 0; j < 8; j++)
            s += Gd[i][j] * q0[j];
        F[i] = s;
    }
}

int AnomalousField::readSigma(string fileWithSigma)
{
    int n;
    ifstream file(fileWithSigma);
    if (!file.is_open())
        return -1;
    file >> n;
    sigmas = new double*[n];
    rotations = new double*[n];
    int num;
    double val;
    for (int i = 0; i < n; i++)
    {
        file >> num;
        sigmas[i] = new double[9];
        for (int j = 0; j < 9; j++)
        {
            file >> val;
            sigmas[num - 1][j] = val;
        }
        rotations[i] = new double[9];
        for (int j = 0; j < 9; j++)
        {
            file >> val;
            rotations[num - 1][j] = val;
        }
    }
    file.close();

    for (int i = 0; i < n; i++)
    {
        double tmp[9]; //MT σ

        // σ = MT σ М
        for (size_t k = 0; k < 3; k++)
        {
            tmp[3 * k] = rotations[i][k] * sigmas[i][0];
            tmp[3 * k + 1] = rotations[i][3 + k] * sigmas[i][4];
            tmp[3 * k + 2] = rotations[i][6 + k] * sigmas[i][8];
        }
        for (size_t k = 0; k < 3; k++)
            for (size_t j = 0; j < 3; j++)
            {
                double sum = 0.0;
                for (size_t l = 0; l < 3; l++)
                    sum += tmp[3 * k + l] * rotations[i][j + 3 * l];
                sigmas[i][3 * k + j] = sum;
            }

        // M = MT M
        for (size_t k = 0; k < 3; k++)
            for (size_t j = 0; j < 3; j++)
            {
                double sum = 0.0;
                for (size_t l = 0; l < 3; l++)
                    sum += rotations[i][k + 3 * l] * rotations[i][j + 3 * l];
                tmp[3 * k + j] = sum;
            }
        for (size_t k = 0; k < 9; k++)
            rotations[i][k] = tmp[k];

    }

    return 0;
}

void AnomalousField::createGlobalSLAE()
{
    double *G0[8][8]; // Локальная матрица жесткости (без сигм)
    double G[8][8]; // Локальная матрица жесткости
    double F[8]; // Локальный вектор правой части

    for (size_t i = 0; i < 8; i++)
        for (size_t j = i; j < 8; j++)
            G0[i][j] = new double[9];

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

        double *sigma = sigmas[nvkat[k]];

        // Создаем локальную матрицу жесткости
        createLocalG(p, G0);

        for (size_t i = 0; i < 8; i++)
            for (size_t j = i; j < 8; j++)
            {
                G[i][j] = 0.0;
                for (size_t l = 0; l < 9; l++)
                    G[i][j] += G0[i][j][l] * sigma[l];
            }


//        for (size_t i = 0; i < 8; i++)
//        {
//            cout << i << ": " << endl;
//            for (size_t j = i; j < 8; j++)
//            {
//                cout << " " << j << ": " << endl;
//                for (size_t l = 0; l < 9; l++)
//                    cout << "  " << G0[i][j][l] << " ";
//                cout << endl;
//            }
//        }


        // Создаем локальный вектор правой части
        createLocalF(p, nvkat[k], G0, F);

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

//    ofstream file("fff");
//    for (size_t i = 0; i < a.n; i++)
//        file << f[i] << endl;
//    file.close();

/*
    f[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1] = 1.0;
    cout << "Координаты источника: "
         << grid->xyz[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1].x << " "
         << grid->xyz[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1].y << " "
         << grid->xyz[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1].z << endl;
*/




    for (size_t i = 0; i < 8; i++)
        for (size_t j = i; j < 8; j++)
            delete[] G0[i][j];
}

void AnomalousField::firstBoundaryCondition()
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
}

void AnomalousField::solve(string method, size_t maxIter, double x0)
{
    // Выделяем память под результат
    if (v)
        delete[] v;
    v = new double[matrix->matrix().n];

    // Начальное приближение
    for (size_t i = 0; i < matrix->matrix().n; i++)
        v[i] = x0;

    // Решаем СЛАУ
    if (method == "MSG_LLT")
        SLAE::solveMSG_LLT(matrix->matrix(), f, v, eps, maxIter);
    else if (method == "LOS_LU")
        SLAE::solveLOS_LU(matrix->matrix(), f, v, eps, maxIter);
}

















