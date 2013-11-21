#include <iostream>
#include "grid.h"

using namespace std;

int main()
{
    Grid2D grid;
    grid.readArea("../normal_field/area.txt");
    grid.readPartitions("../normal_field/grid.txt");
    grid.createGrid();
    grid.saveGrid();
    grid.txtToDat("../normal_field/txtToDat/txttodat");

    return 0;
}

