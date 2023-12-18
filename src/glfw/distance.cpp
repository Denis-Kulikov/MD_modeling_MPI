#include "../include/distance.hpp"

int CompareParticleDistances(const void* a, const void* b) {
    float diff = ((distance_by_index*)b)->dist - ((distance_by_index*)a)->dist;
    return (diff > 0) ? 1 : ((diff < 0) ? -1 : 0);
}

int compare_distances(const void *a, const void *b) {
    double dist_a = *((double *)a);
    double dist_b = *((double *)b);

    if (dist_a < dist_b) {
        return -1;
    } else if (dist_a > dist_b) {
        return 1;
    } else {
        return 0;
    }
}

void sort_distances(const struct distance_by_index &data, int size) {
    // Создаем временный массив для сортировки dist
    double *temp_dist = (double*)malloc(size * sizeof(double));
    if (temp_dist == NULL) {
        // Обработка ошибки выделения памяти
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Копируем значения dist во временный массив
    for (int i = 0; i < size; i++) {
        temp_dist[i] = data.dist[i];
    }

    // Сортируем временный массив
    qsort(temp_dist, size, sizeof(double), compare_distances);

    // Переупорядочиваем index в соответствии с отсортированным dist
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (temp_dist[i] == data.dist[j]) {
                // Обмениваем значения index
                int temp_index = data.index[i];
                data.index[i] = data.index[j];
                data.index[j] = temp_index;
                break;
            }
        }
    }

    // Освобождаем временный массив
    free(temp_dist);
}
