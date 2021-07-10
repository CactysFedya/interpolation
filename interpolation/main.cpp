// interpolation.cpp 
//

#include <iostream>
#include <vector>
#include "interpolation.h"

int main() {
    setlocale(LC_ALL, "Russian");

    std::vector <double> Ox = {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5};
    std::vector <double> Oy = {-2, -1, 0, 0.5, 1, 2.25, 3.5};

    Interpolation object(Ox, Oy);

    object.linear();

    object.canonical();

    object.formulaLagrange();

    object.formulaNewton();
}
