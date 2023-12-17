# MD_modeling_MPI
Molecular dynamics modeling using MPI.

Error: Shader: 0:1(10): error: GLSL 3.30 is not supported. Supported versions are: 1.10, 1.20, 1.30, 1.40, 1.00 ES, and 3.00 ES
Solution: export MESA_GL_VERSION_OVERRIDE=3.3

Error: fatal error: GL/glew.h: No such file or directory #include <GL/glew.h>
Solution: 
sudo apt-get update
sudo apt-get install libglfw3 libglfw3-dev
sudo apt-get install libglew-dev

