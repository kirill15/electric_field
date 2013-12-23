#include <iostream>
#include "normalfield.h"
#include "grid2D.h"
#include <iomanip>

using namespace std;

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

int main()
{
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
    //u0.setSource(1.0);
    cout << "Источник задан" << endl;

    u0.createGlobalRightPart();
    cout << "Глобальный вектор правой части создан" << endl;

    u0.firstBoundaryCondition();
    cout << "Первые краевые учтены" << endl;

    u0.getMatrix()->saveElements();
    cout << "Элементы сохранены" << endl;

    u0.solve();
    cout << "СЛАУ решена" << endl;

    u0.saveSolve();
    cout << "Решене сохранено" << endl;


    unsigned size;
    double *v = u0.getV(size);
    test(v, u0.getGrid()->rz, size);
    cout << "Сравнение с точным напечатано" << endl;


    if (u0.getGrid()->txtToDat("../normal_field/txtToDat/txttodat") == 0)
    {
        cout << "Сетка сохранена в бинарном виде" << endl;

        cout << "Open solve in Telma? (y/n): ";
        char c;
        cin >> c;
        if (c == 'y' or c == 'Y')
            system("cd for_telma && wineconsole go.bat");
    }



/*
    size_t size;
    double *v = u0.getV(size);
    cout << "\nРешение:" << endl;
    for (size_t i = 0; i < size; i++)
        cout << v[i] << endl;
*/


    return 0;
}


/*
Matrix a;
Matrix l;

a.n = 4;

a.ig = new size_t[5];
a.ig[0] = 0; a.ig[1] = 0; a.ig[2] = 1; a.ig[3] = 2; a.ig[4] = 4;

a.jg = new size_t[4];
a.jg[0] = 0; a.jg[1] = 1; a.jg[2] = 1; a.jg[3] = 2;

a.di = new double[4];
a.di[0] = 1; a.di[1] = 101; a.di[2] = 2, a.di[3] = 101;

a.ggl = new double[4];
a.ggl[0] = -1; a.ggl[1] = -10; a.ggl[2] = -10; a.ggl[3] = 1;


SLAE::factorizeLLT(a, l);

cout << "ig: "  << l.ig[0]  << "; " <<  l.ig[1]  << "; " <<  l.ig[2]  << "; " <<  l.ig[3]  << "; " << l.ig[4] << endl;
cout << "jg: "  << l.jg[0]  << "; " <<  l.jg[1]  << "; " <<  l.jg[2]  << "; " <<  l.jg[3]  << endl;
cout << "di: "  << l.di[0]  << "; " <<  l.di[1]  << "; " <<  l.di[2]  << "; " <<  l.di[3]  << endl;
cout << "ggl: " << l.ggl[0] << "; " <<  l.ggl[1] << "; " <<  l.ggl[2] << "; " <<  l.ggl[3] << endl;

return 0;
*/


