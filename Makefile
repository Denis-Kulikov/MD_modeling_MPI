APP_NAME = OpenGL
LIB_NAME = Lib

CC = g++
CFLAGS = -c -Wall -Wextra -Werror
CPPFLAGS = -I src -MP -MMD
LDLIBS = -lglfw -lGL -lGLEW -lm


BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src

APP_PATH = $(BIN_DIR)/$(APP_NAME)

SRC_EXT = cpp
APP_RUN = $(BIN_DIR)/./$(APP_NAME)

APP_SOURCES = $(shell find $(SRC_DIR)/$(APP_NAME) -name '*.$(SRC_EXT)')
APP_OBJECTS = $(APP_SOURCES:$(SRC_DIR)/%.$(SRC_EXT)=$(OBJ_DIR)/%.o)

LIB_SOURCES = $(shell find $(SRC_DIR)/$(LIB_NAME) -name '*.$(SRC_EXT)')
LIB_OBJECTS = $(LIB_SOURCES:$(SRC_DIR)/%.$(SRC_EXT)=$(OBJ_DIR)/%.o)

DEPS = $(APP_OBJECTS:.o=.d) $(LIB_OBJECTS:.o=.d)

.PHONY: all
all: $(APP_PATH)

$(APP_PATH): src/OpenGL/main.cpp src/OpenGL/math_3d.cpp src/OpenGL/pipeline.cpp src/OpenGL/distance.cpp src/OpenGL/n_body.cpp
	$(CC) -o $@ -Wall $^ $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(APP_PATH) $(TEST_PATH) $(LIB_PATH)
	rm -rf $(DEPS) $(APP_OBJECTS) $(LIB_OBJECTS)
	
.PHONY: run
run: $(APP_RUN)
	$(APP_RUN)
