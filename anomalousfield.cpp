#include "anomalousfield.h"

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
double AnomalousField::GetUg(Coord3D p, int nvk)
{
    return 1.0 / (2.0 * M_PI * sqrt(p.x * p.x + p.y *p.y + p.z * p.z) * 0.01);

    //return p.x * p.x * p.y *p.y * p.z * p.z;
    //return exp(p.x * p.y * p.z);
    //return p.x * p.y * p.z;
}

double AnomalousField::getV0(Coord3D p)
{
    double r = sqrt((p.x - pSourcePlus.x) * (p.x - pSourcePlus.x)  +  (p.y - pSourcePlus.y) * (p.y - pSourcePlus.y));
    Coord rz(r, p.z);
    double v0A = v0->getValue(rz);

    return v0A;

//    rz.r = sqrt((p.x - pSourceMinus.x) * (p.x - pSourceMinus.x)  +  (p.y - pSourceMinus.y) * (p.y - pSourceMinus.y));
//    double v0B = v0->getValue(rz);

//    return v0A - v0B;
}

double AnomalousField::getF(Coord3D p, int nvk)
{
    //return -2.0 * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    //return -2.0 * exp(p.x * p.y * p.z) * (p.y*p.y * p.z*p.z  +  p.x*p.x * p.z*p.z  +  p.x*p.x * p.y*p.y);
    return 0.0;
}

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

    Coord3D *xyz = grid->xyz;
    uint *nvkat = grid->nvkat;

    for (size_t i = startIndex; i < limit; i += stepOutside)
        for (size_t j = 0; j < sizeInside; j += stepInside)
        {
            size_t el = i + j;
//            cout << el + 1 << endl;

            if (xyz[el] == Coord3D(0.0, 0.0, 0.0)) continue;

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
    cout << "." << endl;
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

void AnomalousField::setSource(Coord3D plus, Coord3D minus)
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

void AnomalousField::createLocalG(const Coord3D p[8], int nvk, double G[8][8])
{

    double hX = p[1].x - p[0].x;
    double hY = p[2].y - p[0].y;
    double hZ = p[4].z - p[0].z;

//////////////////////////////////////////
/*
    double **sigma = sigmas[nvk];

    double s11hxhy_hz = sigma[0][0] * hX * hY / hZ;
    double s22hxhz_hy = sigma[1][1] * hX * hZ / hY;
    double s33hxhy_hz = sigma[2][2] * hX * hY / hZ;

    G[0][0] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[0][1] = (1.0 / 18.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;
    G[0][2] = (-1.0 / 9.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;
    G[0][3] = (-1.0 / 18.0) * s22hxhz_hy  +  (1.0 / 36.0) * s11hxhy_hz  +  (1.0 / 36.0) * s33hxhy_hz;
    G[0][4] = (1.0 / 18.0) * s22hxhz_hy   -  (1.0 / 9.0) * s11hxhy_hz   -  (1.0 / 9.0) * s33hxhy_hz;
    G[0][5] = (1.0 / 36.0) * s22hxhz_hy   -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;
    G[0][6] = (-1.0 / 18.0) * s22hxhz_hy  -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;
    G[0][7] = (-1.0 / 36.0) * s22hxhz_hy  -  (1.0 / 36.0) * s11hxhy_hz  -  (1.0 / 36.0) * s33hxhy_hz;

    G[1][1] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[1][2] = (-1.0 / 18.0) * s22hxhz_hy  +  (1.0 / 36.0) * s11hxhy_hz  +  (1.0 / 36.0) * s33hxhy_hz;
    G[1][3] = (-1.0 / 9.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;
    G[1][4] = (1.0 / 36.0) * s22hxhz_hy   -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;
    G[1][5] = (1.0 / 18.0) * s22hxhz_hy   -  (1.0 / 9.0) * s11hxhy_hz   -  (1.0 / 9.0) * s33hxhy_hz;
    G[1][6] = (-1.0 / 36.0) * s22hxhz_hy  -  (1.0 / 36.0) * s11hxhy_hz  -  (1.0 / 36.0) * s33hxhy_hz;
    G[1][7] = (-1.0 / 18.0) * s22hxhz_hy  -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;

    G[2][2] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[2][3] = (1.0 / 18.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;
    G[2][4] = (-1.0 / 18.0) * s22hxhz_hy  -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;
    G[2][5] = (-1.0 / 36.0) * s22hxhz_hy  -  (1.0 / 36.0) * s11hxhy_hz  -  (1.0 / 36.0) * s33hxhy_hz;
    G[2][6] = (1.0 / 18.0) * s22hxhz_hy   -  (1.0 / 36.0) * s11hxhy_hz  -  (1.0 / 9.0) * s33hxhy_hz;
    G[2][7] = (1.0 / 36.0) * s22hxhz_hy   -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;

    G[3][3] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[3][4] = (-1.0 / 36.0) * s22hxhz_hy  -  (1.0 / 36.0) * s11hxhy_hz  -  (1.0 / 36.0) * s33hxhy_hz;
    G[3][5] = (-1.0 / 18.0) * s22hxhz_hy  -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;
    G[3][6] = (1.0 / 36.0) * s22hxhz_hy   -  (1.0 / 18.0) * s11hxhy_hz  -  (1.0 / 18.0) * s33hxhy_hz;
    G[3][7] = (1.0 / 18.0) * s22hxhz_hy   -  (1.0 / 9.0) * s11hxhy_hz   -  (1.0 / 9.0) * s33hxhy_hz;

    G[4][4] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[4][5] = (1.0 / 18.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;
    G[4][6] = (-1.0 / 9.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;
    G[4][7] = (-1.0 / 18.0) * s22hxhz_hy  +  (1.0 / 36.0) * s11hxhy_hz  +  (1.0 / 36.0) * s33hxhy_hz;

    G[5][5] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[5][6] = (-1.0 / 18.0) * s22hxhz_hy  +  (1.0 / 36.0) * s11hxhy_hz  +  (1.0 / 36.0) * s33hxhy_hz;
    G[5][7] = (-1.0 / 9.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;

    G[6][6] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
    G[6][7] = (1.0 / 18.0) * s22hxhz_hy   +  (1.0 / 18.0) * s11hxhy_hz  +  (1.0 / 18.0) * s33hxhy_hz;

    G[7][7] = (1.0 / 9.0) * s22hxhz_hy    +  (1.0 / 9.0) * s11hxhy_hz   +  (1.0 / 9.0) * s33hxhy_hz;
*/
//////////////////////////////////////////

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

            G[i][j] =         (Gx[mui][muj] * My[nui][nuj] * Mz[vi][vj] +
                               Mx[mui][muj] * Gy[nui][nuj] * Mz[vi][vj] +
                               Mx[mui][muj] * My[nui][nuj] * Gz[vi][vj]);
        }
    }


}

void AnomalousField::createLocalF(const Coord3D p[], int nvk, double G[][8], double F[])
{
//    double hX = p[1].x - p[0].x;
//    double hY = p[2].y - p[0].y;
//    double hZ = p[4].z - p[0].z;

    double sigma = sigmas[nvk];
    double sigma0 = v0->getSigma((p[4].z + p[0].z) / 2.0);
    double dSigma = sigma0 - sigma;

    // f
    double q0[8];
    for (size_t i = 0; i < 8; i++)
        q0[i] = getV0(p[i]);

    // F=G*q0
    for (size_t i = 0; i < 8; i++)
    {
        double s = 0.0;
        for (size_t j = 0; j < 8; j++)
            s += dSigma * G[i][j] * q0[j];
        F[i] = s;
    }

//    for (size_t i = 0; i < 8; i++)
//        cout << F[i] << endl;
//    cin >> sigma;
}

int AnomalousField::readSigma(string fileWithSigma)
{
    ifstream file(fileWithSigma.c_str());
    if (!file.is_open())
    {
        cerr << "Error with opening " + fileWithSigma;
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


/*
    sigmas = new double**[n];

    for (int i = 0; i < n; i++)
    {
        file >> num;

        sigmas[i] = new double*[3];
        for (size_t j = 0; j < 3; j++)
        {
            sigmas[i] = new double[3];
            for (size_t k = 0; k < 3; k++)
            {
                file >> val;
                sigmas[num - 1][j][k] = val;
            }
        }
    }
*/

    file.close();

    return 0;
}

void AnomalousField::createGlobalSLAE()
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

        double sigma = sigmas[nvkat[k]];

        // Создаем локальную матрицу жесткости
        createLocalG(p, nvkat[k], G);

        // Создаем локальный вектор правой части
///        createLocalF(p, nvkat[k], G, F);

        // Заносим матрицу жесткости в глобальную СЛАУ
        for(int i = 0; i < 8; i++)	// Цикл по строкам локальной матрицы
        {
            a.di[nvtr[k][i]] += sigma * G[i][i];    // Диагональные элементы

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
                a.ggl[left] += sigma * G[i][j];
            }
        }

        // Заносим локальный вектор правой части в глобальный
///        for (size_t i = 0; i < 8; i++)
///            f[nvtr[k][i]] += F[i];
    }

//    cout << "f:\n";
//    for (size_t i = 0; i < a.n; i++)
//        cout << f[i] << endl;

    f[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1] = 1.0 / (2.0 * M_PI);
    cout << "Координаты источника: "
         << grid->xyz[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1].x << " "
         << grid->xyz[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1].y << " "
         << grid->xyz[a.n - (grid->getSizeX() * grid->getSizeY()) / 2 - 1].z << endl;

}

void AnomalousField::firstBoundaryCondition()
{
    Matrix a = matrix->matrix();

    size_t sizeX = grid->getSizeX();
    size_t sizeY = grid->getSizeY();

    size_t sizeXY = sizeX * sizeY;

    // Ближняя грань
    firstBoundaryConditionOnFace(0, sizeXY, sizeX, a.n - sizeXY);
    // Дальняя грань
    firstBoundaryConditionOnFace(sizeXY - sizeX, sizeXY, sizeX, a.n - sizeXY);
    // Левая грань
    firstBoundaryConditionOnFace(sizeX, sizeXY, sizeXY - 2 * sizeX, a.n - sizeXY, sizeX);
    // Правая грань
    firstBoundaryConditionOnFace(2 * sizeX - 1, sizeXY, sizeXY - 2 * sizeX, a.n - sizeXY, sizeX);
    // Нижняя грань
    firstBoundaryConditionOnFace(sizeX + 1, sizeX, sizeX - 2, sizeXY - sizeX - 1);
    // Верхняя грань
    ///firstBoundaryConditionOnFace(a.n - sizeXY + sizeX + 1, sizeX, sizeX - 2, a.n - sizeX);


/*
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
*/
}

void AnomalousField::solve(string method, size_t maxIter)
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
}

















