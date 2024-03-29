#include "../include/glfw.hpp"
#include "../include/MD_modeling.hpp"
#include <sys/time.h> 
#include <omp.h>

#define TRY(command, str_error) try {                                   \
    if (command)                                                        \
        throw std::runtime_error(str_error);                            \
} catch (const std::exception& e) {                                     \
    std::cerr << e.what() << std::endl;                                 \
    abort();                                                            \
}

int nMol;
double size;
DataMol Mol[2];
extern Pipeline pipeline;
Vector3f region;

distance_by_index distances[2];
int iter[2] = {0, 1};
bool IsEnd = false;
bool NextIsEnd = false;

omp_lock_t lock_iter;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void CalculateDistance(distance_by_index &distances, DataMol &Mol)
{
    DO_MOL(nMol) { 
        distances.index[i] = i;
        distances.dist[i] = Mol.p[i].Distance(pipeline.camera.Params.WorldPos);
    }
    sort_distances(distances, nMol);
}

void CloseFiles(FILE **data, int commsize)
{
    for (int i = 0; i < commsize; i++) {
        fclose(data[i]);
    }
    free(data);
}

void ReadData(FILE **data, int commsize, DataMol &Mol)
{
    int NBody = 0;
    int n;

    for (int i = 0; i < commsize; i++) {
        fread(&n, sizeof(int), 1, data[i]);
        fread(&Mol.m[NBody], sizeof(double), n, data[i]);
        fread(&Mol.p[NBody], sizeof(Vector3f), n, data[i]);
        NBody += n;
        IsEnd = feof(data[i]);
    }
}

FILE **ReadParams(int &commsize)
{
    FILE *Params = fopen("data/params.bin", "rb");
    fread(&commsize, sizeof(int), 1, Params);
    fread(&nMol, sizeof(int), 1, Params);
    fread(&size, sizeof(double), 1, Params);
    fclose(Params);

    FILE **data = (FILE**)malloc(sizeof(FILE*) * commsize);
    for (int i = 0; i < commsize; i++) {
        string file_name = "data/result" + std::to_string(i) + ".bin";
        data[i] = fopen(file_name.c_str(), "rb");
    }

    return data;
}

int main(int argc, char** argv)
{
    double t = wtime();
    int commsize;
    FILE **data = ReadParams(commsize);
    region.VSet(size, size, size);
    GLFWwindow* window = nullptr;
    InitializeGLFW(window);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    CompileShaders();
    printf("nMol: %d\n", nMol);

    TRY(((Mol[0].p = (Vector3f*)malloc(sizeof(Vector3f) * nMol)) == nullptr), "Memory allocation error (nMol.p).");
    TRY(((Mol[1].p = (Vector3f*)malloc(sizeof(Vector3f) * nMol)) == nullptr), "Memory allocation error (nMol.p).");
    TRY(((Mol[0].m = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (nMol.m).");
    TRY(((Mol[1].m = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (nMol.m).");
    TRY(((distances[0].index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances[1].index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances[0].dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");
    TRY(((distances[1].dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");

    #pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() == 0) {
            #pragma omp task
            {
                int part = 0;
                int iterCount = 1;
                while (iter[0] == 0);
                while (!IsEnd) {
                    // if ((iterCount % 20) == 0) printf("[%d] %d:%d | %d | DIS=%f\n", iterCount, iter[0], iter[1],  omp_get_thread_num(), distances[part].dist[0]);
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                    RenderSceneCB(distances[part], Mol[part]);
                    glfwSwapBuffers(window);
                    glfwPollEvents();
                    if (iter[part ^ 1] == (iterCount + 1)) {
                        iter[part]++;
                        iterCount++;
                        part ^= 1;
                    }
                }
            }
        } else {
            #pragma omp task
            {
                int part = 0;
                int iterCount = 0;
                while (!IsEnd) {
                    if (iter[part] == iterCount) {
                        ReadData(data, commsize, Mol[part]);
                        ReadData(data, commsize, Mol[part]);
                        CalculateDistance(distances[part], Mol[part]);
                        // if ((iterCount % 20) == 0) printf("[%d] %d:%d | %d\n", iterCount, iter[0], iter[1], omp_get_thread_num());
                        iter[part]++;
                        iterCount++;
                        part ^= 1;
                    }
                }
            }
        }
    }

    CloseFiles(data, commsize);
    free(distances[0].index);
    free(distances[1].index);
    free(distances[0].dist);
    free(distances[1].dist);
    free(Mol[0].m);
    free(Mol[1].m);
    free(Mol[0].p);
    free(Mol[1].p);
    glfwTerminate();

    printf("Total time: %f\n", wtime() - t);

    return 0;
}
