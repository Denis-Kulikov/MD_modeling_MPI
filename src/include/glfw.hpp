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
#include "../include/math_3d.h"
#include "../include/MD_modeling.hpp"
#include "../include/pipeline.hpp"

using namespace std;

#define PI = 3.14159265359f;


void RenderSceneCB(const distance_by_index &distances, const DataMol &Mol);
void CompileShaders();
void InitializeGLFW(GLFWwindow* &window);
