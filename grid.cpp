#include "grid.h"

void Grid2D::createArrays(double *r, double *z, long sizeR, long sizeZ)
{
    // Количество узлов
    countPoints = sizeR * sizeZ;

    // Создаем массив rz
    rz = new Coord[countPoints];
    long k = 0;
    for (long i = 0; i < sizeZ; i++)
    {
        for (long j = 0; j < sizeR; j++)
        {
            rz[k].r = r[j];
            rz[k++].z = z[i];
        }
    }

    // Количество КЭ
    countFE = countPoints - sizeR - sizeZ + 1;

    // Создаем массив nvtr
    nvtr = new long*[countFE];
    k = 0;
    for (long iZ = 0; iZ < sizeZ - 1; iZ++)
    {
        for (long iR = 0; iR < sizeR - 1; iR++)
        {
            nvtr[k] = new long[4];

            nvtr[k][0] = iZ * sizeR + iR;
            nvtr[k][1] = nvtr[k][0] + 1;
            nvtr[k][2] = nvtr[k][0] + sizeR;
            nvtr[k][3] = nvtr[k][2] + 1;

            k++;
        }
    }

    // Создаем массив nvkat
    nvkat = new long[countFE];
    k = 0;
    for (size_t iZ = 0; iZ < partitionsZ.size(); iZ++)                      // Проход по всем областям по оси Z
    {
        for (int i = 0; i < partitionsZ[iZ].first; i++)                     // Проход по всем интервалам по оси Z
        {
            for (size_t iR = 0; iR < partitionsR.size(); iR++)              // Проход по всем областям по оси R
            {
                for (int j = 0; j < partitionsR[iR].first; j++)             // Проход по всем интервалам по оси R
                {
                    nvkat[k++] = areas[partitionsR.size() * iZ + iR][0];
                }
            }
        }
    }
}


Grid2D::Grid2D()
{
}


int Grid2D::readArea(string fileWithArea)
{
    ifstream f(fileWithArea.c_str());
    if (!f.is_open())
    {
        cerr << "File with area is not open";
        return -1;
    }

    f >> nRW;
    rw = new double[nRW];
    for (int i = 0; i < nRW; i++)
        f >> rw[i];

    f >> nZW;
    zw = new double[nZW];
    for (int i = 0; i < nZW; i++)
        f >> zw[i];

    f >> l;
    areas = new int*[l];
    for (int i = 0; i < l; i++)
    {
        areas[i] = new int[5];
        for (int j = 0; j < 5; j++)
        {
            double tmp;
            f >> tmp;
            areas[i][j] = tmp - 1;
        }
    }

    f.close();

    return 0;
}


int Grid2D::readPartitions(string fileWithGrid)
{
    int numberOfIntervals; // Число подыитервалов
    double coefOfDisc; // Коэффициент разрядки

    ifstream f(fileWithGrid.c_str());
    if (!f.is_open())
    {
        cerr << "File with grid is not open";
        return -1;
    }

    // Считываем разбиение каждой подобласти
    for (int i = nRW; i > 1; i--)
    {
        f >> numberOfIntervals >> coefOfDisc;
        if (coefOfDisc < 0)
            coefOfDisc = -1.0 / coefOfDisc;
        partitionsR.push_back(pair<int, double>(numberOfIntervals, coefOfDisc));
    }
    for (int i = nZW; i > 1; i--)
    {
        f >> numberOfIntervals >> coefOfDisc;
        partitionsZ.push_back(pair<int, double>(numberOfIntervals, coefOfDisc));
    }

    f.close();

    return 0;
}


void Grid2D::createGrid()
{
    long sizeR = 1;
    long sizeZ = 1;

    // Размер массива r
    for (size_t i = 0; i < partitionsR.size(); i++)
        sizeR += partitionsR[i].first;

    // Размер массива z
    for (size_t i = 0; i < partitionsZ.size(); i++)
        sizeZ += partitionsZ[i].first;

    // Массивы с координатами
    double *r = new double[sizeR];
    double *z = new double[sizeZ];

    // Заполняем массив r
    int i = 0;
    for (size_t partNum = 0; partNum < partitionsR.size(); partNum++)
    {
        // Левая граница
        double rLeft = rw[areas[partNum][1]];

        // Правая граница
        double rRight = rw[areas[partNum][2]];

        // Количество подобластей
        int numR = partitionsR[partNum].first;

        // Кожффициент разрядки
        double koefRazrR = partitionsR[partNum].second;

        // Начальный шаг
        double stepR = (koefRazrR != 1) ? (rRight - rLeft) * (koefRazrR - 1.0) / (pow(koefRazrR, numR) - 1.0) : (rRight - rLeft) / numR;

        r[i++] = rLeft;
        for (int j = 1; j <= numR; j++, i++)
        {
            r[i] = r[i - 1] + stepR;
            stepR *= koefRazrR;
        }
        r[--i] = rRight;
    }

    // Заполняем массив z
    i = 0;
    for (size_t partNum = 0; partNum < partitionsZ.size(); partNum++)
    {
        // Левая граница
        size_t tmp = partitionsR.size() * partNum;
        double zLow = zw[areas[tmp + partNum][3]];

        // Правая граница
        double zTop = zw[areas[tmp + partNum][4]];

        // Количество подобластей
        int numZ = partitionsZ[partNum].first;

        // Кожффициент разрядки
        double koefRazrZ = partitionsZ[partNum].second;

        // Начальный шаг
        double stepZ = koefRazrZ != 1.0 ? (zTop - zLow) * (koefRazrZ - 1) / (pow(koefRazrZ, numZ) - 1) : (zTop - zLow) / numZ;

        z[i++] = zLow;
        for (int j = 1; j <= numZ; j++, i++)
        {
            z[i] = z[i - 1] + stepZ;
            stepZ *= koefRazrZ;
        }
        z[--i] = zTop;
    }

    // Создаем массивы nvtr, nvkat, rz
    createArrays(r, z, sizeR, sizeZ);

    // Удаляем временные массивы
    delete[] r;
    delete[] z;
}


void Grid2D::saveGrid()
{
    ofstream f;

    f.open("nvtr.txt");
    f << countFE << endl;
    for (long i = 0; i < countFE; i++)
        f << nvtr[i][0] + 1 << " " << nvtr[i][1] + 1 << " " << nvtr[i][2] + 1 << " " << nvtr[i][3] + 1 << endl;
    f.close();

    f.open("nvkat.txt");
    for (long i = 0; i < countFE; i++)
        f << nvkat[i] + 1 << endl;
    f.close();

    f.open("rz.txt");
    f << countPoints << endl;
    for (long i = 0; i < countPoints; i++)
        f << rz[i].r << " " << rz[i].z << endl;
    f.close();
}


int Grid2D::txtToDat(string pathToProgram)
{
    if (execl(pathToProgram.c_str(), "txttodat", NULL) == -1)
    {
        cerr << "Ошибка вызова txttodat" << endl;
        return -1;
    }
    return 0;
}














