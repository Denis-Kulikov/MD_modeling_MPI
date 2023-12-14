cmake_minimum_required(VERSION 3.18.4)

project(
  full-text-search
  VERSION 1.0
  LANGUAGES CXX
)

include_directories(include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

add_subdirectory(src)

add_custom_target(clear
    COMMAND find . -name "CMakeFiles" -exec rm -rf {} +
    COMMAND find . -name "cmake_install.cmake" -exec rm -rf {} +
    COMMAND find . -name "CMakeCache.txt" -exec rm -rf {} +
    COMMAND find . -name "CPackConfig.cmake" -exec rm -rf {} +
    COMMAND find . -name "CPackSourceConfig.cmake" -exec rm -rf {} +
    COMMAND find . -name "DartConfiguration.tcl" -exec rm -rf {} +
    COMMAND find . -name "Makefile" -exec rm -rf {} +
    COMMAND find . -name "bin/MD" -exec rm -rf {} +
    COMMAND rm -rf build
    COMMAND rm -rf build.ninja
    COMMAND rm -rf .ninja_log
    COMMAND rm -rf .ninja_deps
    COMMAND rm -rf _deps
    DEPENDS clean
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Clear"
)
