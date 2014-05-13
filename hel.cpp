#include "hel.h"

void Hel::findNormalField(string fileWithArea, string fileWithGrid, string fileWithSigma)
{
    vZero.setEps(epsVZero);

    cout << "Поиск основного поля:" << endl;

    if (vZero.createGrid(fileWithArea, fileWithGrid) == 0)
        cout << "--Сетка создана" << endl;

    vZero.getGrid()->saveGrid();
    cout << "----Сетка сохранена" << endl;

    vZero.createPortrait();
    cout << "--Портрет сгенерирован" << endl;

    vZero.getMatrix()->savePortrait();
    cout << "----Портрет сохранен" << endl;

    if (vZero.readSigma(fileWithSigma) == 0)
        cout << "--Сигмы считаны" << endl;

    vZero.createGlobalMatrix();
    cout << "--Глобальная СЛАУ собрана" << endl;

    vZero.setSource(J);
    cout << "--Источник задан" << endl;

    vZero.createGlobalRightPart();
    cout << "--Глобальный вектор правой части создан" << endl;

    vZero.firstBoundaryCondition();
    cout << "--Первые краевые учтены" << endl;

    vZero.solve();
    cout << "--СЛАУ решена" << endl;

    vZero.saveSolve();
    cout << "----Решене сохранено" << endl;

    cout << vZero.getValue(Coord(22500.347219543, -200)) << endl;

/*
    if (vZero.getGrid()->txtToDat("../normal_field/txtToDat/txttodat") == 0)
    {
        cout << "Сетка сохранена в бинарном виде" << endl;

        cout << "Open solve in Telma? (y/n): ";
        char c;
        cin >> c;
        if (c == 'y' || c == 'Y')
            system("cd for_telma && wineconsole go.bat");
    }*/
}

void Hel::findAnomalousField(string fileWithArea, string fileWithGrid, string fileWithSigma)
{
    vPlus.setEps(epsVPlus);

    if (vPlus.createGrid(fileWithArea, fileWithGrid, 2) == 0)
        cout << "--Сетка создана" << endl;

    vPlus.getGrid()->saveGrid();
    cout << "----Сетка сохранена" << endl;

    vPlus.createPortrait();
    cout << "--Портрет сгенерирован" << endl;

    vPlus.getMatrix()->savePortrait();
    cout << "----Портрет сохранен" << endl;

    vPlus.setNormalField(&vZero);

    vPlus.setSource(anode, cathode);

    if (vPlus.readSigma(fileWithSigma) == 0)
        cout << "--Сигмы считаны" << endl;

    vPlus.createGlobalSLAE();
    cout << "--Глобальная СЛАУ собрана" << endl;

    vPlus.firstBoundaryCondition();
    cout << "--Первые краевые учтены" << endl;

    vPlus.getMatrix()->saveElements();

    vPlus.solve();
    cout << "--СЛАУ решена" << endl;


    ofstream file("uu.csv");
    size_t size;
    double *x = vPlus.getV(size);
    Coord3D *p = vPlus.getGrid()->xyz;
    double sum = 0.0, sumX = 0.0;
    for (size_t i = 0; i < size; i++)
    {
        if (p[i].x == 0 && p[i].y == 0) continue;
        double uS = 1.0 / (2.0 * M_PI * sqrt(p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z) * 0.01);
        file << p[i].x << ";" << p[i].y << ";" << p[i].z << ";" << x[i] << ";" << uS << ";" << x[i] - uS << endl;
        sum += (x[i] - uS) * (x[i] - uS);
        sumX += uS * uS;
    }

    file << "||u*-u|| / ||u*||:;" << sqrt(sum) / sqrt(sumX) << endl;

    file.close();

    cout << "( 25000.0,  25000.0, 0.0)\t" << vPlus.getValue(Coord3D(25000.0, 25000.0, 0.0)) << endl;
    cout << "( 25000.0, -25000.0, 0.0)\t" << vPlus.getValue(Coord3D(25000.0, -25000.0, 0.0)) << endl;
    cout << "(-25000.0, -25000.0, 0.0)\t" << vPlus.getValue(Coord3D(-25000.0, -25000.0, 0.0)) << endl;
    cout << "(-25000.0,  25000.0, 0.0)\t" << vPlus.getValue(Coord3D(-25000.0, 25000.0, 0.0)) << endl;
    cout << "(22500.0, -125.0, -200.0)\t" << vPlus.getValue(Coord3D(22500.0, -125.0, -200.0)) << endl;
}

void Hel::setJ(double J)
{
    this->J = J;
}

void Hel::setEps(double epsForV0, double epsForVPlus)
{
    epsVZero = epsForV0;
    epsVPlus = epsForVPlus;
}

void Hel::setHelCoords(Coord3D anode, Coord3D cathode)
{
    this->anode = anode;
    this->cathode = cathode;
}

Hel::Hel()
{
}
