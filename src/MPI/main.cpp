#include "../include/glfw.hpp"
#include <mpi.h>

extern Pipeline pipeline;
extern struct distance_by_index distances;
extern Vector3f region;
extern Vector3i initUcell;
extern DataMol Mol;
extern int nMol, moreCycles, stepCount, stepLimit;
extern int width, height;

#define FIRST(x, s) (x * nMol / s)
#define LAST(x, s) (FIRST(x, s) + nMol / s)
#define NELEMS(x, s) (LAST(x, s) - FIRST(x, s) + (nMol % s) * (x == (s - 1)))

MPI_Datatype vector3f_type;

void Alltoall(int rank, int size)
{ 
    // printf("[%d] Alltoall\n", rank);
    MPI_Request reqs[(size - 1) * 2 * 4];
    MPI_Status stats[(size - 1) * 2 * 4];
    int req_count = 0;

    for (int i = 0; i < size; i++)
    {
        if (i != rank) {
            // MPI_Isend(Mol.p + FIRST(i, size), NELEMS(i, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            // MPI_Irecv(Mol.p + FIRST(rank, size), NELEMS(rank, size), vector3f_type, rank, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Isend(Mol.p + FIRST(rank, size), NELEMS(rank, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(Mol.p + FIRST(i, size), NELEMS(i, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);

            MPI_Isend(Mol.v + FIRST(rank, size), NELEMS(rank, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(Mol.v + FIRST(i, size), NELEMS(i, size), vector3f_type, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);

            MPI_Isend(distances.dist + FIRST(rank, size), NELEMS(rank, size), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(distances.dist + FIRST(i, size), NELEMS(i, size), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);

            MPI_Isend(distances.index + FIRST(rank, size), NELEMS(rank, size), MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(distances.index + FIRST(i, size), NELEMS(i, size), MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[req_count++]);


            // printf("[%d] Isend: begin add(%d) size(%d)\tIrecv: begin add(%d) size(%d)\n", i,
            // FIRST(i, size),
            // NELEMS(i, size),
            // FIRST(rank, size), 
            // NELEMS(rank, size));
        } 
    }

    MPI_Waitall((size - 1) * 2, reqs, stats);
}
    // for (int i = 0; i < size - 1; i++) {
    //     for (int j = i + 1; j < size; j++) {
    //         MPI_Isend(Mol.p + FIRST(i, size), LAST(i, size) - FIRST(i, size) + (nMol % size) * (i == (size - 1)), vector3f_type, j, 0, MPI_COMM_WORLD, &reqs[req_count++]);
    //         MPI_Irecv(Mol.p + FIRST(j, size), LAST(j, size) - FIRST(j, size) + (nMol % size) * (j == (size - 1)), vector3f_type, j, 0, MPI_COMM_WORLD, &reqs[req_count++]);
    //     }
    // }
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
    GLFWwindow* window = nullptr;
    MPI_Init(&argc, &argv);
    double ttotal = -MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm cartcomm;

    MPI_Type_contiguous(3, MPI_DOUBLE, &vector3f_type);
    MPI_Type_commit(&vector3f_type);

    SetParams();
    TRY(!AllocArrays(), "Memory allocation error (AllocArrays).");
    TRY(((distances.index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances.dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");

    // for (int i = 0; i < nMol; i++) {
    //     distances->dist = i;
    //     distances->index = i;
    // }

    if (rank == 0) {
        InitializeGLFW(window);
        glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
        CompileShaders();
        
        SetupJob ();
        printf("Body: %d\n\n", nMol);
        Vector3f CameraPos(0.0f, 0.1f, -region.z * 2 - pow(2, region.z));
        Vector3f CameraTarget(0.0f, 0.0f, 1.0f);
        Vector3f CameraUp(0.0f, 1.0f, 0.0f);
        pipeline.camera.SetCamera(CameraPos, CameraTarget, CameraUp);
        pipeline.camera.SetPerspectiveProj(60.0f, width, height, 0.5f, 1000.0f);
        pipeline.object.SetScale(0, 0, 0);
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
        Alltoall(rank, commsize);
        if (rank == 0) {
            // printf("Last: %f %f %f\n", Mol.p[last + 1].x, Mol.p[last + 1].y, Mol.p[last + 1].z);
        // Передача между процессами 

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            RenderSceneCB();
            glfwSwapBuffers(window);
            glfwPollEvents();
            if (stepCount >= stepLimit) moreCycles = 0;
        }
    }

    // if (rank == 1) {
        // printf("RC: %d %d\n", first, last);
        // for (int i = 0; i < 4; i++)
            // printf("coord: %f %f %f\n", Mol.p[i].x, Mol.p[i].y, Mol.p[i].z);
    // }

    ttotal += MPI_Wtime();
    printf("ttotal: %f\n", ttotal);

    // MPI_Type_free(&vector3f_type);

    // free(distances);
    // free(Mol.m);
    // free(Mol.v);
    // free(Mol.f);
    // free(Mol.p);

    if (rank == 0) {
        glfwTerminate();
    }

    MPI_Finalize();

    return 0;
}
    // printf ("%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
    //         stepCount, timeNow, VCSum (vSum) / nMol, PropEst (totEnergy), 5
    //         PropEst (kinEnergy), PropEst (pressure));
