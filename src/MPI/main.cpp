// #include "../include/glfw.hpp"
#include "../include/MD_modeling.hpp"
#include "../include/math_3d.h"
#include <mpi.h>

#define FIRST(x, s)  (static_cast<int>(nMol / s) * x + std::min(x, nMol % s))
#define LAST(x, s)   (FIRST(x, s) + static_cast<int>(nMol / s) + (x < nMol % s ? 1 : 0))
#define NELEMS(x, s) (LAST(x, s) - FIRST(x, s) + (nMol % s) * (x == (s - 1)))

extern DataMol Mol;
extern Vector3f vSum, Total_vSum;
extern double uSum, virSum, vvSum, Total_uSum, Total_virSum, Total_vvSum;
extern int nMol, moreCycles, stepCount, stepLimit, stepAvg;

MPI_Datatype vector3f_type;

void ExchangeAndReduce(int rank, int size)
{
    MPI_Request reqs[(size - 1) * 2];
    MPI_Status stats[(size - 1) * 2];
    int req_count = 0;

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Isend(Mol.p + FIRST(rank, size), NELEMS(rank, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]); // передача позиций точек
            MPI_Irecv(Mol.p + FIRST(i, size), NELEMS(i, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
        } 
    }

    if (((stepCount + 1) % stepAvg) == 0) {
        Total_uSum = 0, Total_virSum = 0, Total_vvSum = 0;
        Total_vSum.VZero();

        MPI_Reduce(&Total_uSum, &uSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&Total_virSum, &virSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&Total_vvSum, &vvSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&Total_vSum.x, &vSum.x, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    MPI_Waitall((size - 1) * 2, reqs, stats);
}

void BroadcastDataMol(int rank, int size)
{
    MPI_Bcast(Mol.p, nMol, vector3f_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(Mol.v, nMol, vector3f_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(Mol.m, nMol, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char** argv)
{
    int commsize, rank;
    MPI_Init(&argc, &argv);
    double ttotal = -MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Type_contiguous(3, MPI_DOUBLE, &vector3f_type);
    MPI_Type_commit(&vector3f_type);

    SetParams();
    TRY(!AllocArrays(), "Memory allocation error (AllocArrays).");

    if (rank == 0) {
        SetupJob ();
        printf("Body: %d\n", nMol);
        printf("stepLimit: %d\n\n", stepLimit);
    }
    BroadcastDataMol(rank, commsize);
    MPI_Barrier(MPI_COMM_WORLD);

    int first, last;
    first = FIRST(rank, commsize);
    last = LAST(rank, commsize);
    printf("[%d]: %d %d\n", rank, first, last);

    while (moreCycles) {
        SingleStep (first, last);
        ExchangeAndReduce(rank, commsize);
        if (stepCount >= stepLimit) moreCycles = 0;
    }

    ttotal += MPI_Wtime();
    printf("[%d] ttotal: %f\n", rank, ttotal);

    MPI_Type_free(&vector3f_type);

    free(Mol.m);
    free(Mol.v);
    free(Mol.f);
    free(Mol.p);

    MPI_Finalize();

    return 0;
}
