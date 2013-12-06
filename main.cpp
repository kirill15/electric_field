#include <iostream>
#include "normalfield.h"
#include "grid2D.h"

using namespace std;

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

    u0.createGlobalMatrix();
    cout << "Глобальная СЛАУ собрана" << endl;

    u0.setSource(1, 1);
    cout << "Источник задан" << endl;

    u0.createGlobalRightPart();;
    cout << "Глобальный вектор правой части создан" << endl;

    u0.getMatrix()->saveElements();
    cout << "Элементы сохранены" << endl;

    u0.solve();
    cout << "СЛАУ решена" << endl;

    ofstream f;
    f.open("v.txt");
    size_t size;
    double *v = u0.getV(size);
    for (int i = 0; i < size; ++i) {
        f << v[i] << endl;
    }
    f.close();

    if (u0.getGrid()->txtToDat("../normal_field/txtToDat/txttodat") == 0)
        cout << "Сетка сохранена в бинарном виде" << endl;


/*
    size_t size;
    double *v = u0.getV(size);
    cout << "\nРешение:" << endl;
    for (size_t i = 0; i < size; i++)
        cout << v[i] << endl;
*/

    return 0;
}

