#include "../include/glfw.hpp"

extern int nMol;
extern double size;
extern DataMol Mol;
extern FILE *result;
extern struct distance_by_index distances;
extern Pipeline pipeline;
extern Vector3f region;

bool IsEnd = false;

void CalculateDistance()
{
    DO_MOL(nMol) { 
        distances.index[i] = i;
        distances.dist[i] = Mol.p[i].Distance(pipeline.camera.Params.WorldPos);
    }
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
    result = fopen("data/result.bin", "rb");
    ReadParams();
    region.VSet(size, size, size);
    GLFWwindow* window = nullptr;
    InitializeGLFW(window);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    CompileShaders();

    TRY(((Mol.p = (Vector3f*)malloc(sizeof(Vector3f) * nMol)) == nullptr), "Memory allocation error (nMol.p).");
    TRY(((Mol.m = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (nMol.m).");
    TRY(((distances.index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances.dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");

    ReadMass();
    ReadPosition();

    while (!IsEnd) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        CalculateDistance();
        RenderSceneCB();
        ReadPosition();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    fclose(result);
    free(distances.index);
    free(distances.dist);
    free(Mol.m);
    free(Mol.p);
    glfwTerminate();

    return 0;
}
