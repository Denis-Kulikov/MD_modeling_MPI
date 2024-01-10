APP_NAME = MD
DEBUG_NAME = ${APP_NAME}_debug
MPI_NAME = ${APP_NAME}_MPI

CC = g++
MPICXX = mpicxx
CFLAGS = -c -Wall -Wextra -Werror -MMD -MP
LDLIBS = -lglfw -lGL -lGLEW -lm


BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src

APP_PATH = $(BIN_DIR)/$(APP_NAME)
DEBUG_PATH = $(BIN_DIR)/$(DEBUG_NAME)
MPI_PATH = $(BIN_DIR)/$(MPI_NAME)

SRC_EXT = cpp
APP_RUN = $(BIN_DIR)/./$(APP_NAME)
DEBUG_RUN = $(BIN_DIR)/./$(DEBUG_NAME)
MPI_RUN = $(BIN_DIR)/./$(MPI_NAME)

APP_SOURCES = $(shell find $(SRC_DIR)/$(APP_NAME) -name '*.$(SRC_EXT)')
APP_OBJECTS = $(APP_SOURCES:$(SRC_DIR)/%.$(SRC_EXT)=$(OBJ_DIR)/%.o)

LIB_SOURCES = $(shell find $(SRC_DIR)/$(LIB_NAME) -name '*.$(SRC_EXT)')
LIB_OBJECTS = $(LIB_SOURCES:$(SRC_DIR)/%.$(SRC_EXT)=$(OBJ_DIR)/%.o)

DEPS = $(APP_OBJECTS:.o=.d) $(LIB_OBJECTS:.o=.d)

.PHONY: all
all: $(APP_PATH)

$(APP_PATH): src/glfw/main.cpp src/glfw/glfw.cpp src/Math/math_3d.cpp src/glfw/pipeline.cpp src/glfw/distance.cpp
	$(CC)  -fopenmp -o $@ -Wall $^ $(LDLIBS)

.PHONY: debug
debug: $(DEBUG_PATH)

$(DEBUG_PATH): src/glfw/main.cpp src/glfw/glfw.cpp src/Math/math_3d.cpp src/glfw/pipeline.cpp src/glfw/distance.cpp 
	$(CC) -g -fopenmp -o $@ -Wall $^ $(LDLIBS)
	
.PHONY: mpi
mpi: $(MPI_PATH)

$(MPI_PATH): src/MPI/all-cell.cpp  src/Math/math_3d.cpp src/MD_modeling/MD_modeling.cpp
	$(MPICXX) -o $@ -Wall $^ -lm
	
.PHONY: clean
clean:
	rm -f $(APP_PATH) $(TEST_PATH) $(LIB_PATH)
	rm -rf $(DEPS) $(APP_OBJECTS) $(LIB_OBJECTS)
	
.PHONY: run
run: $(APP_RUN)
	$(APP_RUN)
	
.PHONY: mrun
mrun: $(MPI_RUN)
	mpiexec $(MPI_RUN)
	
.PHONY: drun
drun: $(DEBUG_RUN)
	gdb $(DEBUG_RUN)