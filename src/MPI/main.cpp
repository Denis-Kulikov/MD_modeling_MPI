#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <unistd.h>

#include "../include/MD_modeling.hpp"
#include "../include/math_3d.h"

int main(int argc, char** argv)
{
    init();

    while (moreCycles) {
        SingleStep ();
        if (stepCount >= stepLimit) moreCycles = 0;
    }

    // printf ("%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
    //         stepCount, timeNow, VCSum (vSum) / nMol, PropEst (totEnergy), 5
    //         PropEst (kinEnergy), PropEst (pressure));

    free(Mol.m);
    free(Mol.v);
    free(Mol.f);
    free(Mol.p);

    return 0;
}
