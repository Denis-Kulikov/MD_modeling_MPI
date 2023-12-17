APP_NAME = MD
LIB_NAME = Lib

CC = g++
MPICXX = mpicxx
CFLAGS = -c -Wall -Wextra -Werror
CPPFLAGS = -I src -MP -MMD
LDLIBS = -lglfw -lGL -lGLEW -lm


BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src

APP_PATH = $(BIN_DIR)/$(APP_NAME)
DEBUG_PATH = $(BIN_DIR)/$(APP_NAME)_debug
MPI_PATH = $(BIN_DIR)/$(APP_NAME)_mpi

SRC_EXT = cpp
APP_RUN = $(BIN_DIR)/./$(APP_NAME)

APP_SOURCES = $(shell find $(SRC_DIR)/$(APP_NAME) -name '*.$(SRC_EXT)')
APP_OBJECTS = $(APP_SOURCES:$(SRC_DIR)/%.$(SRC_EXT)=$(OBJ_DIR)/%.o)

LIB_SOURCES = $(shell find $(SRC_DIR)/$(LIB_NAME) -name '*.$(SRC_EXT)')
LIB_OBJECTS = $(LIB_SOURCES:$(SRC_DIR)/%.$(SRC_EXT)=$(OBJ_DIR)/%.o)

DEPS = $(APP_OBJECTS:.o=.d) $(LIB_OBJECTS:.o=.d)

.PHONY: all
all: $(APP_PATH)

$(APP_PATH): src/OpenGL/main.cpp src/Math/math_3d.cpp src/OpenGL/pipeline.cpp src/OpenGL/distance.cpp src/MD_modeling/MD_modeling.cpp
	$(CC) -o $@ -Wall $^ $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(APP_PATH) $(TEST_PATH) $(LIB_PATH)
	rm -rf $(DEPS) $(APP_OBJECTS) $(LIB_OBJECTS)
	
.PHONY: run
run: $(APP_RUN)
	$(APP_RUN)

.PHONY: debug
debug: $(DEBUG_PATH)

$(DEBUG_PATH): src/OpenGL/main.cpp src/Math/math_3d.cpp src/OpenGL/pipeline.cpp src/OpenGL/distance.cpp src/MD_modeling/MD_modeling.cpp
	$(CC) -g -o $@ -Wall $^ $(LDLIBS)
	
.PHONY: drun
drun: $(DEBUG_PATH)
	gdb $(DEBUG_PATH)
	
.PHONY: mpi
debug: $(DEBUG_PATH)

$(MPI_PATH): src/OpenGL/main.cpp src/Math/math_3d.cpp src/MD_modeling/MD_modeling.cpp
	$(MPICXX) -g -o $@ -Wall $^ -lm
	
.PHONY: drun
drun: $(DEBUG_PATH)
	gdb $(DEBUG_PATH)