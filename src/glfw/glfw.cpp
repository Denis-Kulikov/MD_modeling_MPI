// #include <fstream>
// #include <iostream>
// #include <cmath>
// #include <cassert>
// #include <sstream>
// #include <string.h>
// #include <stdio.h>
// #include <vector>
// #include <unistd.h>

// #include <mpi.h>

// #include <GL/glew.h>
// #include <GLFW/glfw3.h>

// #include "../include/distance.hpp"
// #include "../include/MD_modeling.hpp"
// #include "../include/math_3d.h"
// #include "../include/pipeline.hpp"
// #include "../include/pipeline.hpp"

#include "../include/glfw.hpp"

GLuint VAO;
GLuint VBO;
GLuint EBO;
GLuint gWorldLocation;
GLuint gScaleLocation;
GLuint gMassLocation;
GLuint ShaderSphere;
GLuint ShaderCube;

int width = 1280;
int height = 768;

double radius = 0.075f;
float PI = 3.14159265359f;

int SPHERE_SEGMENTS = 16;
int numVertices;
int numIndices;

Pipeline pipeline;
struct distance_by_index distances;
extern Vector3f region;
extern DataMol Mol;
extern int nMol, moreCycles, stepCount, stepLimit;

void CreateCube()
{
    std::vector<GLfloat> vertices = {
        -1, -1, 1,
        -1, 1, 1,   -1, 1, -1,  -1, 1, 1, 
        1, 1, 1,    1, 1, -1,    1, 1, 1,
        1, -1, 1,   1, -1, -1,   1, -1, 1,
        -1, -1, 1,

        -1, -1, -1,
        -1, 1, -1,
        1, 1, -1,
        1, -1, -1,
        -1, -1, -1
    };

    numVertices = static_cast<int>(vertices.size()) / 3;

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void createSphere()
{
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;

    for (int lat = 0; lat <= SPHERE_SEGMENTS; ++lat) {
        double theta = lat * M_PI / SPHERE_SEGMENTS;
        double sinTheta = sin(theta);
        double cosTheta = cos(theta);

        for (int lon = 0; lon <= SPHERE_SEGMENTS; ++lon) {
            double phi = lon * 2 * M_PI / SPHERE_SEGMENTS;
            double sinPhi = sin(phi);
            double cosPhi = cos(phi);

            float x = static_cast<float>(cosPhi * sinTheta);
            float y = static_cast<float>(cosTheta);
            float z = static_cast<float>(sinPhi * sinTheta);

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
    }

    numVertices = static_cast<int>(vertices.size()) / 3;

    for (int lat = 0; lat < SPHERE_SEGMENTS; ++lat) {
        for (int lon = 0; lon < SPHERE_SEGMENTS; ++lon) {
            int current = lat * (SPHERE_SEGMENTS + 1) + lon;
            int next = current + SPHERE_SEGMENTS + 1;

            indices.push_back(current);
            indices.push_back(next);
            indices.push_back(current + 1);

            indices.push_back(current + 1);
            indices.push_back(next);
            indices.push_back(next + 1);
        }
    }

    numIndices = static_cast<int>(indices.size());

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void DrawCube()
{
    pipeline.object.SetWorldPos(0, 0, 0);
    pipeline.object.SetRotate(0, 0, 0);
    pipeline.object.SetScale(region.x + radius, region.y + radius, region.z + radius);

    GLfloat floatMatrix[16];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) 
            floatMatrix[i * 4 + j] = static_cast<GLfloat>(pipeline.GetTrans()->m[i][j]);
    }

    glUniformMatrix4fv(gWorldLocation, 1, GL_TRUE, floatMatrix);
    glBindVertexArray(VAO);
    glDrawArrays(GL_LINE_STRIP, 0, numVertices);
    glBindVertexArray(0);
}

void DrawSphere(int i)
{
    pipeline.object.SetWorldPos(Mol.p[i].x, Mol.p[i].y, Mol.p[i].z);

    GLfloat floatMatrix[16];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) 
            floatMatrix[i * 4 + j] = static_cast<GLfloat>(pipeline.GetTrans()->m[i][j]);
    }

    glUniformMatrix4fv(gWorldLocation, 1, GL_TRUE, floatMatrix);
    glUniform1f(gMassLocation, static_cast<float>(Mol.m[i]));

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void RenderSceneCB()
{
    sort_distances(distances, nMol);
    // qsort(distances, nMol, sizeof(*distances), CompareParticleDistances);
    glUseProgram(ShaderSphere); 
    createSphere();
    pipeline.object.SetScale(radius, radius, radius);
    pipeline.object.SetRotate(0, 0, 0);
    for (int i = 0; i < nMol; i++) {
        int particleIndex = distances.index[i];
        DrawSphere(particleIndex);
    }

    glUseProgram(ShaderCube);
    CreateCube();
    DrawCube();
}

GLuint LoadShader(const char *shader_path, GLuint type)
{
    ifstream shader_file(shader_path);

    if (!shader_file.is_open()) {
        cerr << "Error: Could not open shader file '" << shader_path << "'" << endl;
        return 0;
    }

    stringstream shader_stream;
    shader_stream << shader_file.rdbuf();
    shader_file.close();

    char* code = (char*)malloc(shader_stream.str().length() + 1);
    for (size_t i = 0; i < shader_stream.str().length(); ++i) {
        code[i] = shader_stream.str()[i];
    }
    code[shader_stream.str().length()] = '\0';
    const GLchar* shader_code = code;

    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &shader_code, NULL);
    glCompileShader(shader);

    GLint ok;
    GLchar log[2000];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        glGetShaderInfoLog(shader, 2000, NULL, log);
        printf("Shader(%s): %s\n", shader_path, log);
        cout << shader_stream.str().c_str() << endl;
    }

    return shader;
}

void CompileShadersProgram(GLuint ShaderProgram, const char *FS, const char *VS)
{
    GLuint shader_color;
    GLuint shader_position;

    if (FS != nullptr) {
        shader_color = LoadShader(FS, GL_FRAGMENT_SHADER);
        glAttachShader(ShaderProgram, shader_color);
    }
    if (VS != nullptr) {
        shader_position = LoadShader(VS, GL_VERTEX_SHADER);
        glAttachShader(ShaderProgram, shader_position);
    }

    glLinkProgram(ShaderProgram);

    GLint ok;
    GLchar log[2000];
    glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &ok);
    if (!ok) {
        glGetProgramInfoLog(ShaderProgram, 2000, NULL, log);
        std::cout << "ShaderSphere Compilation Log:\n" << log << std::endl;
    
        if (FS != nullptr) {
            GLint infoLogLength;
            glGetShaderiv(shader_color, GL_INFO_LOG_LENGTH, &infoLogLength);
            GLchar *infoLog = new GLchar[infoLogLength + 1];
            glGetShaderInfoLog(shader_color, infoLogLength, NULL, infoLog);
            std::cout << "ShaderSphere shader_color Log:\n" << infoLog << std::endl;
            delete[] infoLog;
        }
        if (VS != nullptr) {
            GLint infoLogLength;
            glGetShaderiv(shader_position, GL_INFO_LOG_LENGTH, &infoLogLength);
            GLchar *infoLog = new GLchar[infoLogLength + 1];
            glGetShaderInfoLog(shader_position, infoLogLength, NULL, infoLog);
            std::cout << "ShaderSphere shader_position Log:\n" << infoLog << std::endl;
            delete[] infoLog;
        }
    }

    glUseProgram(ShaderProgram);
}

void CompileShaders()
{
    ShaderSphere = glCreateProgram();
    ShaderCube   = glCreateProgram();

    CompileShadersProgram(ShaderSphere, "shaders/sphere_fs.glsl", "shaders/sphere_vs.glsl");
    gScaleLocation = glGetUniformLocation(ShaderSphere, "gWorld");
    gMassLocation = glGetUniformLocation(ShaderSphere, "gMass");
    assert((gScaleLocation != 0xFFFFFFFF) || (gMassLocation != 0xFFFFFFFF));

    CompileShadersProgram(ShaderCube, nullptr, "shaders/cube_vs.glsl");
    gScaleLocation = glGetUniformLocation(ShaderCube, "gWorld");
    assert(gScaleLocation != 0xFFFFFFFF);
}

static void KeyboardCB(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    double speed_movement = 0.1 * region.z;
    double speed_rotation = 0.075;
    switch (key) {
        case GLFW_KEY_F:
            glfwSetWindowShouldClose(window, GLFW_TRUE);
            moreCycles = 0;
            break;
        case GLFW_KEY_W:
            pipeline.camera.Params.WorldPos.z += speed_movement;
            break;
        case GLFW_KEY_S:
            pipeline.camera.Params.WorldPos.z -= speed_movement;
            break;
        case GLFW_KEY_D:
            pipeline.camera.Params.WorldPos.x += speed_movement;
            break;
        case GLFW_KEY_A:
            pipeline.camera.Params.WorldPos.x -= speed_movement;
            break;
        case GLFW_KEY_SPACE:
            pipeline.camera.Params.WorldPos.y += speed_movement;
            break;
        case GLFW_KEY_C:
            pipeline.camera.Params.WorldPos.y -= speed_movement;
            break;
        case GLFW_KEY_E:
            pipeline.camera.Params.Target.x += speed_rotation;
            break;
        case GLFW_KEY_Q:
            pipeline.camera.Params.Target.x -= speed_rotation;
            break;
    }
}

void InitializeGLFW(GLFWwindow* &window)
{
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(width, height, "NBody", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, KeyboardCB);

    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }

    glEnable(GL_DEPTH_TEST);

    glClearColor(0.1f, 0.1f, 0.1f, 0.0f);
}
