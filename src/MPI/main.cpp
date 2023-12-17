// #include "../include/glfw.hpp"
#include "../include/math_3d.h"
#include "../include/MD_modeling.hpp"
#include <mpi.h>

// extern struct distance_by_index distances;
extern Vector3f region;
extern Vector3i initUcell;
extern DataMol Mol;
extern int nMol, moreCycles, stepCount, stepLimit;
extern int width, height;
extern double uSum, virSum;

#define FIRST(x, s) (x * nMol / s)
#define LAST(x, s) (FIRST(x, s) + nMol / s)
#define NELEMS(x, s) (LAST(x, s) - FIRST(x, s) + (nMol % s) * (x == (s - 1)))

MPI_Datatype vector3f_type;

void ExchangeAndReduce(int rank, int size)
{
    MPI_Request reqs[(size - 1) * 2 * 2];
    MPI_Status stats[(size - 1) * 2 * 2];
    int req_count = 0;

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Isend(Mol.p + FIRST(rank, size), NELEMS(rank, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(Mol.p + FIRST(i, size), NELEMS(i, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);

            MPI_Isend(Mol.v + FIRST(rank, size), NELEMS(rank, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(Mol.v + FIRST(i, size), NELEMS(i, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);

            // MPI_Isend(distances.dist + FIRST(rank, size), NELEMS(rank, size), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            // MPI_Irecv(distances.dist + FIRST(i, size), NELEMS(i, size), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);

            // MPI_Isend(distances.index + FIRST(rank, size), NELEMS(rank, size), MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            // MPI_Irecv(distances.index + FIRST(i, size), NELEMS(i, size), MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
        } 
    }
    double NewuSum = 0, NewvirSum = 0;
    MPI_Reduce(&NewuSum, &uSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&NewvirSum, &virSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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
    // GLFWwindow* window = nullptr;
    MPI_Init(&argc, &argv);
    double ttotal = -MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Type_contiguous(3, MPI_DOUBLE, &vector3f_type);
    MPI_Type_commit(&vector3f_type);

    SetParams();
    TRY(!AllocArrays(), "Memory allocation error (AllocArrays).");
    // TRY(((distances.index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    // TRY(((distances.dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");

    if (rank == 0) {
        // InitializeGLFW(window);
        // glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
        // CompileShaders();
        SetupJob ();
        printf("Body: %d\n", nMol);
        printf("stepLimit: %d\n\n", stepLimit);
    }
    BroadcastDataMol(rank, commsize);
    MPI_Barrier(MPI_COMM_WORLD);

    int first, last;
    if (rank == (commsize - 1)) {
        first = FIRST(rank, commsize);
        last = LAST(rank, commsize) + nMol % commsize;
    } else {
        first = FIRST(rank, commsize);
        last = LAST(rank, commsize);
    }
    printf("FL: %d %d\n", first, last);

    while (moreCycles) {
        SingleStep (first, last);
        ExchangeAndReduce(rank, commsize);
        if (rank == 0) {
            // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            // RenderSceneCB();
            // glfwSwapBuffers(window);
            // glfwPollEvents();
        }
        if (stepCount >= stepLimit) moreCycles = 0;
    }

    ttotal += MPI_Wtime();
    printf("[%d] ttotal: %f\n", rank, ttotal);

    MPI_Type_free(&vector3f_type);

    // free(distances.index);
    // free(distances.dist);
    free(Mol.m);
    free(Mol.v);
    free(Mol.f);
    free(Mol.p);

    // if (rank == 0) {
        // glfwTerminate();
    // }

    MPI_Finalize();

    return 0;
}
