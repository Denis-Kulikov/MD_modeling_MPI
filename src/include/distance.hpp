#pragma once
#include <stdlib.h>
#include <iostream>

struct distance_by_index {
    int *index;
    double *dist;
};

int CompareParticleDistances(const void* a, const void* b);
void sort_distances(const struct distance_by_index &data, int size);
