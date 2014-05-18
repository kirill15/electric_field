#include "hel.h"

using namespace std;


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

//    cout << vZero.getValue(Coord(22500.347219543, -200)) << endl;

/*
    if (vZero.getGrid()->txtToDat("../normal_field/txtToDat/txttodat") == 0)
    {
        cout << "Сетка сохранена в бинарном виде" << endl;

        cout << "Open solve in Telma? (y/n): ";
        char c;
        cin >> c;
        if (c == 'y' || c == 'Y')
            system("cd for_telma && wineconsole go.bat");
    }
    */


/*
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    for (int i = 0; i < 15000; i+=10)
        for (int j = 0; j < 9999; j+=10)
            cout << vZero.getValue(Coord(i+1.0, -j-1.0)) << endl;
    gettimeofday(&t2, NULL);
    cout << t2.tv_sec - t1.tv_sec << endl;
    cout << t2.tv_usec - t1.tv_usec << endl;
    return;
*/



    cout << "10, 10, -50:\t" << vZero.getValue(Coord(sqrt(10*10 + 10*10), -50)) << endl;
    cout << "10, 10, -150:\t" << vZero.getValue(Coord(sqrt(10*10 + 10*10), -150)) << endl;
    cout << "10, 10, -250:\t" << vZero.getValue(Coord(sqrt(10*10 + 10*10), -250)) << endl;
    cout << "70, 70, -200:\t" << vZero.getValue(Coord(sqrt(70*70 + 70*70), -200)) << endl;
    cout << "200, 200, -150:\t" << vZero.getValue(Coord(sqrt(200*200 + 200*200), -150)) << endl;
    cout << "200, 200, -250:\t" << vZero.getValue(Coord(sqrt(200*200 + 200*200), -250)) << endl;
    cout << "150, 150, -1000:\t" << vZero.getValue(Coord(sqrt(150*150 + 150*150), -1000)) << endl;

    ofstream file("V0.csv");
    size_t size;
    double *x = vZero.getV(size);
    Coord *p = vZero.getGrid()->rz;
    for (size_t i = 0; i < size; i++)
    {
        file<< p[i].r << ";" << p[i].z << ";" << x[i] << endl;
    }
    file.close();
}

void Hel::findAnomalousField(string fileWithArea, string fileWithGrid, string fileWithSigma)
{
    vPlus.setEps(epsVPlus);

    cout << "--Создание сетки" << endl;
    if (vPlus.createGrid(fileWithArea, fileWithGrid, 0) == 0)
        cout << "----ОК (" << vPlus.getGrid()->getCountFE() << " К.Э.)" << endl;

//    vPlus.getGrid()->saveGrid();
//    cout << "----Сетка сохранена" << endl;

    cout << "--Генерация портрета" << endl;
    vPlus.createPortrait();
    cout << "----OK" << endl;

//    vPlus.getMatrix()->savePortrait();
//    cout << "----Портрет сохранен" << endl;

    vPlus.setNormalField(&vZero);

    vPlus.setSource(anode, cathode);

    if (vPlus.readSigma(fileWithSigma) != 0)
        cout << "--Ошибка чтения сигм!" << endl;

    cout << "--Сборка глобальной СЛАУ" << endl;
    vPlus.createGlobalSLAE();
    cout << "----OK" << endl;

//    vPlus.getMatrix()->saveElements();

    cout << "--Учет краевых условий" << endl;
    vPlus.firstBoundaryCondition();
    cout << "----OK" << endl;

//    vPlus.getMatrix()->saveElements();

    cout << "--Решение СЛАУ" << endl;
//    vPlus.solve("LOS_LU", 2000);
    vPlus.solve("MSG_LLT", 2000);
    cout << "----OK" << endl;








/*
    cout << "10, 10, -50:\t" << vPlus.getValue(Coord3D(10, 10, -50)) << endl;
    cout << "10, 10, -150:\t" << vPlus.getValue(Coord3D(10, 10, -150)) << endl;
    cout << "10, 10, -250:\t" << vPlus.getValue(Coord3D(10, 10, -250)) << endl;
    cout << "70, 70, -200:\t" << vPlus.getValue(Coord3D(70, 70, -200)) << endl;
    cout << "200, 200, -150:\t" << vPlus.getValue(Coord3D(200, 200, -150)) << endl;
    cout << "200, 200, -250:\t" << vPlus.getValue(Coord3D(200, 200, -250)) << endl;
    cout << "150, 150, -1000:\t" << vPlus.getValue(Coord3D(150, 150, -1000)) << endl;
*/





/*
    ofstream file("V+.csv");
    size_t size;
    double *x = vPlus.getV(size);
    Coord3D *p = vPlus.getGrid()->xyz;
//    double sum = 0.0, sumX = 0.0;
    for (size_t i = 0; i < size; i++)
    {
//        if (p[i] == Coord3D(0, 0, 0)) continue;
//        double uS = 1.0 / (2.0 * M_PI * sqrt(p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z) * 0.01);
        file << p[i].x << ";" << p[i].y << ";" << p[i].z << ";" << x[i] << endl;//<< ";" << uS << ";" << x[i] - uS << endl;
//        sum += (x[i] - uS) * (x[i] - uS);
//        sumX += uS * uS;
    }

//    file << "||u*-u|| / ||u*||:;" << sqrt(sum) / sqrt(sumX) << endl;

    file.close();
*/

/*
    double xx[5];
    xx[0] = fabs(vPlus.getValue(Coord3D(100.0, 100.0, 0.0)) - 1.0 / (2.0 * M_PI * sqrt(100*100 + 100*100) * 0.01));
    cout << "( 100.0, 100.0,  0.0)\t" << vPlus.getValue(Coord3D(100.0, 100.0, 0.0)) << "\t" << 1.0 / (2.0 * M_PI * sqrt(100*100 + 100*100) * 0.01) << "  \t " << std::scientific << xx[0] << std::fixed << endl;
    xx[1] = fabs(vPlus.getValue(Coord3D(200.0, 200.0, 0.0)) - 1.0 / (2.0 * M_PI * sqrt(200*200 + 200*200) * 0.01));
    cout << "( 200.0, 200,0,  0.0)\t" << vPlus.getValue(Coord3D(200.0, 200.0, 0.0)) << "\t" << 1.0 / (2.0 * M_PI * sqrt(200*200 + 200*200) * 0.01) << "\t " << std::scientific << xx[1] << std::fixed << endl;
    xx[2] = fabs(vPlus.getValue(Coord3D(300.0, 300.0, 0.0)) - 1.0 / (2.0 * M_PI * sqrt(300*300 + 300*300) * 0.01));
    cout << "(300.0,  300.0,  0.0)\t" << vPlus.getValue(Coord3D(300.0, 300.0, 0.0)) << "\t" << 1.0 / (2.0 * M_PI * sqrt(300*300 + 300*300) * 0.01) << "\t " << std::scientific << xx[2] << std::fixed << endl;
    xx[3] = fabs(vPlus.getValue(Coord3D(500.0, 500.0, 0.0)) - 1.0 / (2.0 * M_PI * sqrt(500*500 + 500*500) * 0.01));
    cout << "(500.0,  500.0,  0.0)\t" << vPlus.getValue(Coord3D(500.0, 500.0, 0.0)) << "\t" << 1.0 / (2.0 * M_PI * sqrt(500*500 + 500*500) * 0.01) << "\t " << std::scientific << xx[3] << std::fixed << endl;
    xx[4] = fabs(vPlus.getValue(Coord3D(1000.0, 1000.0, 0.0)) - 1.0 / (2.0 * M_PI * sqrt(1000*1000 + 1000*1000) * 0.01));
    cout << "(1000.0, 1000.0, 0.0)\t" << vPlus.getValue(Coord3D(1000.0, 1000.0, 0.0)) << "\t" << 1.0 / (2.0 * M_PI * sqrt(1000*1000 + 1000*1000) * 0.01) << "\t " << std::scientific << xx[4] << std::fixed << endl;
    cout << "(-8189.87, 6068.66, -200)\t" << vPlus.getValue(Coord3D(-8189.87, 6068.66, -200)) << "\t" << 1.0 / (2.0 * M_PI * sqrt(8189.87*8189.87 + 6068.66*6068.66 + 200*200) * 0.01) << endl;

    cout << "Норма: " << std::scientific << sqrt(xx[0]*xx[0] + xx[1]*xx[1] + xx[2]*xx[2] + xx[3]*xx[3] + xx[4]*xx[4]) << endl;

    cout << "(10, 10, -147.419): " << vPlus.getValue(Coord3D(10, 10, -147.419)) << endl;*/
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

void Hel::setHelCoords(const Coord3D &anode, const Coord3D &cathode)
{
    this->anode = anode;
    this->cathode = cathode;
}

double Hel::getValue(Coord3D p, bool isNotAnomalous)
{
    return vZero.getValue(Coord(sqrt((p.x - anode.x) * (p.x - anode.x)  +  (p.y - anode.y) * (p.y - anode.y)), p.z))
         - vZero.getValue(Coord(sqrt((p.x - cathode.x) * (p.x - cathode.x)  +  (p.y - cathode.y) * (p.y - cathode.y)), p.z))
         + ((isNotAnomalous) ? 0.0 : vPlus.getValue(p));
}

NormalField *Hel::getNormalField()
{
    return &vZero;
}

AnomalousField *Hel::getAnomalousField()
{
    return &vPlus;
}

Hel::Hel()
{
}







