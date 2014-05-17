#ifndef HEL_H
#define HEL_H

#include "normalfield.h"
#include "anomalousfield.h"

//#include "times.h"
#include <sys/time.h>

class Hel
{
private:
    NormalField vZero; // Основное поле
    AnomalousField vPlus; // Добавочное поле

    double J; // Плотность тока источника (модуль)

    double epsVZero; // Точность решения СЛАУ для v0
    double epsVPlus; // Точность решения СЛАУ для v+

    Coord3D anode;
    Coord3D cathode;


public:
    Hel();

    // Найти основное поле
    void findNormalField(string fileWithArea, string fileWithGrid, string fileWithSigma);

    // Найти добавочное поле
    void findAnomalousField(string fileWithArea, string fileWithGrid, string fileWithSigma);

    // Задать плотность источника
    void setJ(double J);

    // Задать точности решения СЛАУ
    void setEps(double epsForV0, double epsForVPlus);

    // Задать координаты анода и катода
    void setHelCoords(const Coord3D &anode, const Coord3D &cathode);

    // Значение V в точке
    double getValue(Coord3D p);

    NormalField *getNormalField();

    AnomalousField *getAnomalousField();
};

#endif // HEL_H
