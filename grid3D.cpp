#include "grid3D.h"

unsigned Grid3D::getAreaNumber(int xLeft, int xRight, int yNear, int yFar, int zLow, int zTop)
{
    for (uint i = 0; i < l; i++)
        if (areas[i][1] <= xLeft && areas[i][2] >= xRight && areas[i][3] <= yNear && areas[i][4] >= yFar && areas[i][5] <= zLow && areas[i][6] >= zTop)
            return areas[i][0];
    throw "Not find number of area";
}


void Grid3D::createArrays(double *x, double *y, double *z)
{
    // Количество узлов
    countPoints = sizeX * sizeY * sizeZ;

    // Создаем массив xyz
    xyz = new Coord3D[countPoints];
    unsigned k = 0;
    for (size_t i = 0; i < sizeZ; i++)
        for (size_t j = 0; j < sizeY; j++)
            for (size_t g = 0; g < sizeX; g++)
            {
                xyz[k].x = x[g];
                xyz[k].y = y[j];
                xyz[k++].z = z[i];
            }

    // Количество КЭ
    countFE = (sizeX - 1) * (sizeY - 1) * (sizeZ - 1);

    // Создаем массив nvtr
    nvtr = new uint*[countFE];
    k = 0;
    for (size_t iZ = 0; iZ < sizeZ - 1; iZ++)
        for (size_t iY = 0; iY < sizeY - 1; iY++)
            for (size_t iX = 0; iX < sizeX - 1; iX++)
            {
                nvtr[k] = new uint[8];

                nvtr[k][0] = iZ * sizeX * sizeY + iY * sizeX + iX;
                nvtr[k][1] = nvtr[k][0] + 1;
                nvtr[k][2] = nvtr[k][0] + sizeX;
                nvtr[k][3] = nvtr[k][2] + 1;

                nvtr[k][4] = nvtr[k][0] + sizeX * sizeY;
                nvtr[k][5] = nvtr[k][4] + 1;
                nvtr[k][6] = nvtr[k][4] + sizeX;
                nvtr[k][7] = nvtr[k][6] + 1;

                k++;
            }

    // Создаем массив nvkat
    nvkat = new uint[countFE];
    k = 0;
    for (size_t iZ = 0; iZ < partitionsZ.size(); iZ++)                      // Проход по всем областям по оси Z
        for (int i = 0; i < partitionsZ[iZ].first; i++)                     // Проход по всем интервалам по оси Z
            for (size_t iY = 0; iY < partitionsY.size(); iY++)              // Проход по всем областям по оси Y
                for (int j = 0; j < partitionsY[iY].first; j++)             // Проход по всем интервалам по оси Y
                    for (size_t iX = 0; iX < partitionsX.size(); iX++)      // Проход по всем областям по оси X
                    {
                        uint num = getAreaNumber(iX, iX + 1, iY, iY + 1, iZ, iZ + 1);
                        for (int g = 0; g < partitionsZ[iX].first; g++)     // Проход по всем интервалам по оси X
                            nvkat[k++] = num;
                    }
}


Grid3D::Grid3D()
{
}


int Grid3D::readArea(string fileWithArea)
{
    ifstream f(fileWithArea.c_str());
    if (!f.is_open())
    {
        cerr << "File with area is not open";
        return -1;
    }

    f >> nXW;
    xw = new double[nXW];
    for (int i = 0; i < nXW; i++)
        f >> xw[i];

    f >> nYW;
    yw = new double[nYW];
    for (int i = 0; i < nYW; i++)
        f >> yw[i];

    f >> nZW;
    zw = new double[nZW];
    for (int i = 0; i < nZW; i++)
        f >> zw[i];

    f >> l;
    areas = new int*[l];
    for (size_t i = 0; i < l; i++)
    {
        areas[i] = new int[7];
        for (size_t j = 0; j < 7; j++)
        {
            double tmp;
            f >> tmp;
            areas[i][j] = tmp - 1;
        }
    }

    f.close();

    return 0;
}


int Grid3D::readPartitions(string fileWithGrid, size_t fragmentation)
{
    int numberOfIntervals; // Число подъитервалов
    double coefOfDisc; // Коэффициент разрядки

    ifstream f(fileWithGrid.c_str());
    if (!f.is_open())
    {
        cerr << "File with grid is not open";
        return -1;
    }

    // Считываем разбиение каждой подобласти
    for (int i = nXW; i > 1; i--)
    {
        f >> numberOfIntervals >> coefOfDisc;
        if (coefOfDisc < 0)
            coefOfDisc = -1.0 / coefOfDisc;

        // Дробим сетку fragmentation раз
        for (size_t j = 0; j < fragmentation; j++)
        {
            numberOfIntervals *= 2;
            coefOfDisc = sqrt(coefOfDisc);
        }

        partitionsX.push_back(pair<int, double>(numberOfIntervals, coefOfDisc));
    }

    for (int i = nYW; i > 1; i--)
    {
        f >> numberOfIntervals >> coefOfDisc;
        if (coefOfDisc < 0)
            coefOfDisc = -1.0 / coefOfDisc;

        // Дробим сетку fragmentation раз
        for (size_t j = 0; j < fragmentation; j++)
        {
            numberOfIntervals *= 2;
            coefOfDisc = sqrt(coefOfDisc);
        }

        partitionsY.push_back(pair<int, double>(numberOfIntervals, coefOfDisc));
    }

    for (int i = nZW; i > 1; i--)
    {
        f >> numberOfIntervals >> coefOfDisc;
        if (coefOfDisc < 0)
            coefOfDisc = -1.0 / coefOfDisc;

        // Дробим сетку fragmentation раз
        for (size_t j = 0; j < fragmentation; j++)
        {
            numberOfIntervals *= 2;
            coefOfDisc = sqrt(coefOfDisc);
        }

        partitionsZ.push_back(pair<int, double>(numberOfIntervals, coefOfDisc));
    }

    f.close();

    return 0;
}


void Grid3D::createGrid()
{
    sizeX = sizeY = sizeZ = 1;

    // Размер массива x
    for (size_t i = 0; i < partitionsX.size(); i++)
        sizeX += partitionsX[i].first;

    // Размер массива y
    for (size_t i = 0; i < partitionsY.size(); i++)
        sizeY += partitionsY[i].first;

    // Размер массива z
    for (size_t i = 0; i < partitionsZ.size(); i++)
        sizeZ += partitionsZ[i].first;

    // Массивы с координатами
    double *x = new double[sizeX];
    double *y = new double[sizeY];
    double *z = new double[sizeZ];

    // Заполняем массив x
    int i = 0;
    size_t numInArea = 0;
    for (size_t partNum = 0; partNum < partitionsX.size(); partNum++)
    {
        // Левая граница
        double xLeft = xw[areas[numInArea][1]];

        // Правая граница
        double xRight = xw[areas[numInArea++][2]];

        // Количество подобластей
        int numX = partitionsX[partNum].first;

        // Коэффициент разрядки
        double koefRazrX = partitionsX[partNum].second;

        // Начальный шаг
        double stepX = (koefRazrX != 1) ? (xRight - xLeft) * (koefRazrX - 1.0) / (pow(koefRazrX, numX) - 1.0) : (xRight - xLeft) / numX;

        x[i++] = xLeft;
        for (int j = 1; j <= numX; j++, i++)
        {
            x[i] = x[i - 1] + stepX;
            stepX *= koefRazrX;
        }
        x[--i] = xRight;
    }

    // Заполняем массив y
    i = 0;
    for (size_t partNum = 0; partNum < partitionsY.size(); partNum++)
    {
        // Ближняя граница
        double yLeft = yw[areas[numInArea][3]];

        // Дальняя граница
        double yRight = yw[areas[numInArea++][4]];

        // Количество подобластей
        int numY = partitionsY[partNum].first;

        // Коэффициент разрядки
        double koefRazrY = partitionsY[partNum].second;

        // Начальный шаг
        double stepY = (koefRazrY != 1) ? (yRight - yLeft) * (koefRazrY - 1.0) / (pow(koefRazrY, numY) - 1.0) : (yRight - yLeft) / numY;

        y[i++] = yLeft;
        for (int j = 1; j <= numY; j++, i++)
        {
            y[i] = y[i - 1] + stepY;
            stepY *= koefRazrY;
        }
        y[--i] = yRight;
    }

    // Заполняем массив z
    i = 0;
    for (size_t partNum = 0; partNum < partitionsZ.size(); partNum++)
    {
        // Нижняя граница
        double zLow = zw[areas[numInArea][5]];

        // Верхняя граница
        double zTop = zw[areas[numInArea++][6]];

        // Количество подобластей
        int numZ = partitionsZ[partNum].first;

        // Коэффициент разрядки
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

    // Создаем массивы nvtr, nvkat, xyz
    createArrays(x, y, z);

    // Удаляем временные массивы
    delete[] x;
    delete[] y;
    delete[] z;
}
