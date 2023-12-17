#include "../include/glfw.hpp"

extern int nMol;
extern double size;
extern DataMol Mol;
extern FILE *result;
extern struct distance_by_index distances;

bool IsEnd = false;

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
    fread(&size, sizeof(int), 1, result);
}

int main(int argc, char** argv)
{
    GLFWwindow* window = nullptr;
    InitializeGLFW(window);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    CompileShaders();
    result = fopen("data/result.bin", "rb");
    ReadParams();

    cout << "Beg" << endl;

    TRY(((Mol.p = (Vector3f*)malloc(sizeof(Vector3f) * nMol)) == nullptr), "Memory allocation error (nMol.p).");
    TRY(((Mol.m = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (nMol.m).");
    TRY(((distances.index = (int*)malloc(sizeof(int) * nMol)) == nullptr), "Memory allocation error (distances.index).");
    TRY(((distances.dist = (double*)malloc(sizeof(double) * nMol)) == nullptr), "Memory allocation error (distances.dist).");

    cout << "Mem" << endl;
    ReadMass();
    ReadPosition();
    cout << nMol << "\t" << size << endl;

    while (!IsEnd) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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
