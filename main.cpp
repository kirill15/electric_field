#include <iostream>
#include <iomanip>
#include "grid3D.h"
#include "hel.h"

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
    //return p.x * p.x * p.y *p.y * p.z * p.z;
    //return exp(p.x * p.y * p.z);
    return p.x * p.y * p.z;
}



int main()
{
    /*
    NormalField u0;

    if (u0.createGrid("../normal_field/area.txt", "../normal_field/grid.txt") == 0)
        cout << "Сетка создана" << endl;

    u0.getGrid()->saveGrid();
    cout << "Сетка сохранена в текстовом виде" << endl;

    u0.createPortrait();
    cout << "Портрет сгенерирован" << endl;

    u0.getMatrix()->savePortrait();
    cout << "Портрет сохранен" << endl;

    if (u0.readSigma("../normal_field/sigma.txt") == 0)
        cout << "Сигмы считаны" << endl;

    u0.createGlobalMatrix();
    cout << "Глобальная СЛАУ собрана" << endl;

    u0.setSource(1.0 / (2.0 * M_PI));
    cout << "Источник задан" << endl;

    u0.createGlobalRightPart();
    cout << "Глобальный вектор правой части создан" << endl;

    u0.firstBoundaryCondition();
    cout << "Первые краевые учтены" << endl;

    u0.getMatrix()->saveElements();
    cout << "Элементы сохранены" << endl;

    u0.solve("LOS_LU");
    cout << "СЛАУ решена" << endl;

    u0.saveSolve();
    cout << "Решене сохранено" << endl;

    cout << u0.getValue(Coord(10, -10)) << endl;

    cout << u0.getCountFE() << " конечных элементов." << endl;

    if (u0.getGrid()->txtToDat("../normal_field/txtToDat/txttodat") == 0)
    {
        cout << "Сетка сохранена в бинарном виде" << endl;

        cout << "Open solve in Telma? (y/n): ";
        char c;
        cin >> c;
        if (c == 'y' or c == 'Y')
            system("cd for_telma && wineconsole go.bat");
    }
    */


/*
    Mke3D u;
    u.createGrid("../normal_field/area3D.txt", "../normal_field/grid3D.txt", 0);
    u.getGrid()->saveGrid();
    u.createPortrait();
    u.getMatrix()->savePortrait();
    u.readSigma("../normal_field/sigma3D.txt");
    u.createGlobalSLAE();
    //u.getMatrix()->saveElements();
    u.firstBoundaryCondition();
    u.getMatrix()->saveElements();
    u.setEps(1e-15);
    u.solve();




    cout << "( 100.0, 100.0,  -1.0)\t" << u.getValue(Coord3D(100.0, 100.0, -1.0)) << endl;
    cout << "( 200.0, 200,0,  -1.0)\t" << u.getValue(Coord3D(200.0, 200.0, -1.0)) << endl;
    cout << "(300.0,  300.0,  -1.0)\t" << u.getValue(Coord3D(300.0, 300.0, -1.0)) << endl;
    cout << "(500.0,  500.0,  -1.0)\t" << u.getValue(Coord3D(500.0, 500.0, -1.0)) << endl;
    cout << "(1000.0, 1000.0, -1.0)\t" << u.getValue(Coord3D(1000.0, 1000.0, -1.0)) << endl;
    cout << "(-8189.87, 6068.66, -1)\t" << u.getValue(Coord3D(-8189.87, 6068.66, -1)) << endl;



    double *x, *xe;
    size_t size;
    x = u.getV(size);
    Coord3D *pp = u.getGrid()->xyz;
    xe = new double[size];
    cout << "x:\n";
    for (size_t i = 0; i < size; i++)
    {
        xe[i] = getXYZExact(pp[i]);
        cout << i + 1 << ":\t" << x[i] << "\t" << xe[i] << endl;
    }

    double normXE = SLAE::norm(xe, size);
    SLAE::subVectorVector(x, xe, size, xe);
    double normX_Xe = SLAE::norm(xe, size);
    cout << "К.Э: "<< std::scientific << u.getGrid()->getCountFE() << endl;
    cout << "Относительная погрешность: " << normX_Xe / normXE << endl;
*/



    Hel hel;
    hel.setJ(1.0 / (2.0 * M_PI));
    hel.setEps(1e-15, 1e-15);
    hel.setHelCoords(Coord3D(0, 0, 0), Coord3D(0, 0, 0));
    hel.findNormalField("../normal_field/area.txt", "../normal_field/grid.txt", "../normal_field/sigma.txt");
    hel.findAnomalousField("../normal_field/area3D.txt", "../normal_field/grid3D.txt", "../normal_field/sigma3D.txt");




    return 0;
}




/*
    unsigned size;
    double *v = u0.getV(size);
    test(v, u0.getGrid()->rz, size);
    cout << "Сравнение с точным напечатано" << endl;
*/
