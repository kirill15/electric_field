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

    if (u0.getGrid()->txtToDat("../normal_field/txtToDat/txttodat") == 0)
        cout << "Сетка сохранена в бинарном виде" << endl;

    u0.createPortrait();
    cout << "Портрет сгенерирован" << endl;

    u0.getMatrix()->savePortrait();
    cout << "Портрет сохранен" << endl;

    //cout << u0.getGrid()->nvtr[0][0];

    //return 1;

    u0.createGlobalSLAE();
    cout << "Глобальная СЛАУ собрана" << endl;

    u0.getMatrix()->saveElements();
    cout << "Элементы сохранены" << endl;

    return 0;
}

