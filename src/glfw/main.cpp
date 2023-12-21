#include "../include/glfw.hpp"
#include <sys/time.h> 
#include <omp.h>

extern int nMol;
extern double size;
extern DataMol Mol;
extern FILE *result;
extern Pipeline pipeline;
extern Vector3f region;

distance_by_index distances[2];
int iter[2] = {0, 1};
bool IsEnd = false;

omp_lock_t lock_iter;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void CalculateDistance(distance_by_index &distances)
{
    DO_MOL(nMol) { 
        distances.index[i] = i;
        distances.dist[i] = Mol.p[i].Distance(pipeline.camera.Params.WorldPos);
    }
    sort_distances(distances, nMol);
}

void ReadPosition()
{
    DO_MOL(nMol) fread(&Mol.p[i], sizeof(Vector3f), 1, result);
    IsEnd = feof(result);
}

void ReadMass()
{
    DO_MOL(nMol) fread(&Mol.m[i], sizeof(double), 1, result);
}

void ReadParams()
{
    fread(&nMol, sizeof(int), 1, result);
    fread(&size, sizeof(double), 1, result);
}

int main(int argc, char** argv)
{
    double t = wtime();
    TRY(((result = fopen("data/result.bin", "rb")) == nullptr), "No such file: data/result.bin\n");
    ReadParams();
    region.VSet(size, size, size);
    GLFWwindow* window = nullptr;
    InitializeGLFW(window);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    CompileShaders();
    printf("nMol: %d\n", nMol);

    TRY(((Mol.p = (Vector3f*)malloc(sizeof(Vector3f) * nMol)) == nullptr), "Memory allocation error (nMol.p).");
    TRY(((Mol.m = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (nMol.m).");
    TRY(((distances[0].index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances[1].index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances[0].dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");
    TRY(((distances[1].dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");

    ReadMass();

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
                    RenderSceneCB(distances[part]);
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
                        ReadPosition();
                        CalculateDistance(distances[part]);
                        // if ((iterCount % 20) == 0) printf("[%d] %d:%d | %d\n", iterCount, iter[0], iter[1], omp_get_thread_num());
                        iter[part]++;
                        iterCount++;
                        part ^= 1;
                    }
                }
            }
        }
    }

    fclose(result);
    free(distances[0].index);
    free(distances[1].index);
    free(distances[0].dist);
    free(distances[1].dist);
    free(Mol.m);
    free(Mol.p);
    glfwTerminate();

    printf("Total time: %f\n", wtime() - t);

    return 0;
}
