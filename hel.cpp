#include "hel.h"

void hel::findNormalField(string fileWithArea, string fileWithGrid, string fileWithSigma)
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
}

void hel::findAnomalousField(string fileWithArea, string fileWithGrid, string fileWithSigma)
{
    vPlus.setEps(epsVPlus);

    if (vPlus.createGrid(fileWithArea, fileWithGrid) == 0)
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

    vPlus.solve();
    cout << "--СЛАУ решена" << endl;
}

hel::hel()
{
}
