#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <unistd.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "../include/distance.hpp"
#include "../include/MD_modeling.hpp"
#include "../include/math_3d.h"
#include "../include/pipeline.hpp"

using namespace std;

// void CreateCube();
// void createSphere();
// void DrawCube();
// void DrawSphere(int i);
void RenderSceneCB();
// GLuint LoadShader(const char *shader_path, GLuint type);
// void CompileShadersProgram(GLuint ShaderProgram, const char *FS, const char *VS);
void CompileShaders();
// static void KeyboardCB(GLFWwindow* window, int key, int scancode, int action, int mods);
void InitializeGLFW(GLFWwindow* &window);
