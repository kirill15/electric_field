#include <iostream>
#include <iomanip>
#include "hel.h"
#include <unistd.h>

//#include "mke_3d.h"

using namespace std;


// Сравнение решения с точным (для двумерной области)
void test(double *v, Coord *rz, size_t size)
{
    double sigma = 0.01;

    long double pogr = 0.0;
    long double normVStar = 0.0;

    ofstream file("compare.txt");
    file << setw(15) << "r" << setw(15) << "z" << setw(15) << "Vчисл.\t" << setw(15) << "Vточн.\t" << setw(15) << "|Vчисл. - Vточн.|" << endl;
    for (size_t i = 0; i < size; i+=1)
    {
        double v_star = 1.0 / (2.0 * M_PI * sqrt(rz[i].r * rz[i].r + rz[i].z * rz[i].z) * sigma);
        file << setw(15) << rz[i].r << setw(15) << rz[i].z << setw(15) << v[i] << setw(15) << v_star << setw(15) << fabs(v[i] - v_star) << endl;
        pogr += (rz[i].r + rz[i].z == 0) ? 0 : fabs(v[i] - v_star) * fabs(v[i] - v_star);
        if (rz[i].r + rz[i].z != 0)
            normVStar += v_star * v_star;
    }
    file << "Норма вектора погрешностей: " << setw(15) << sqrt(pogr) / sqrt(normVStar);
    file.close();
}

double getXYZExact(Coord3D p)
{
    //return p.x * p.y * p.z;
    //return exp(p.x * p.y * p.z);
    //return p.x*p.x * p.y*p.y * p.z*p.z;
    //return p.x*p.x*p.x * p.y*p.y*p.y * p.z*p.z*p.z;
    return p.x*p.x*p.x*p.x * p.y*p.y*p.y*p.y * p.z*p.z*p.z*p.z;
}



int main()
{
    double opt_epsForNormalField = 1e-15;
    double opt_epsForAnomalousField = 1e-15;
    double opt_J = 1.0 / (2.0 * M_PI);
    double opt_helX = 200.0;

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    Hel hel;
    hel.setJ(opt_J);
    hel.setEps(opt_epsForNormalField, opt_epsForAnomalousField);
    hel.setHelCoords(Coord3D(opt_helX, 0, 0), Coord3D(-opt_helX, 0, 0));
    hel.findNormalField("../normal_field/area.txt", "../normal_field/grid.txt", "../normal_field/sigma.txt");
    hel.findAnomalousField("../normal_field/area3D.txt", "../normal_field/grid3D.txt", "../normal_field/sigma3D.txt");

    gettimeofday(&t2, NULL);
    cout << "Время: " << (long long)t2.tv_sec * 1000 + t2.tv_usec / 1000 - (long long)t1.tv_sec * 1000 + t1.tv_usec / 1000 << " мс" << endl;

    ofstream graf("grafff_1");
    bool isNOTanomalous = false;

//    int koord1 = -8000, koord2 = koord1 + 50;
//    while (koord2 <= 8000)
//    {
//        graf << std::fixed << (koord2 + koord1) / 2.0 << "\t" << std::scientific << hel.getValue(Coord3D(koord1, 0, 0), isNOTanomalous) - hel.getValue(Coord3D(koord2, 0, 0), isNOTanomalous) << endl;

//        koord1 += 100;
//        koord2 = koord1 + 50;
//    }

//    size_t tmp;
//    cout << "0, 0: " << hel.getNormalField()->getValue(Coord(0, 0)) << endl;


//    for (size_t i = 0; i < 100; i++)
//        graf << std::fixed  << 100 * i << "\t" << std::scientific << hel.getNormalField()->getValue(Coord(fabs(100.0 * i - 200.0), 0)) << endl;


    hel.CreateSectionXZ(0.0);
    if (Grid2D::txtToDat("../normal_field/txtToDat/txttodat") == 0)
    {
        cout << "Сетка сохранена в бинарном виде" << endl;

        cout << "Open solve in Telma? (y/n): ";
        char c;
        cin >> c;
        if (c == 'y' || c == 'Y')
            system("cd for_telma && wineconsole go.bat");
    }
    else
        cout << "Err" << endl;

    cout << std::fixed  << 200 << "\t" << std::scientific << hel.getValue(Coord3D(200, 0, 0)) << endl;


    return 0;
}




/*
    unsigned size;
    double *v = u0.getV(size);
    test(v, u0.getGrid()->rz, size);
    cout << "Сравнение с точным напечатано" << endl;
*/
