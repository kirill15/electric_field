#ifndef HEL_H
#define HEL_H

#include "normalfield.h"
#include "anomalousfield.h"

class hel
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

    // Найти основное поле
    void findNormalField(string fileWithArea, string fileWithGrid, string fileWithSigma);

    // Найти добавочное поле
    void findAnomalousField(string fileWithArea, string fileWithGrid, string fileWithSigma);


    hel();
};

#endif // HEL_H
