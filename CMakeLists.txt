cmake_minimum_required(VERSION 3.14)
project(CruiseMissileSim LANGUAGES CXX)

# 默认生成器: MinGW64 (MinGW Makefiles)。见 CMakePresets.json 的 mingw64 preset。
#   命令行构建:  cmake --preset mingw64 && cmake --build build
#   VSCode:      使用 .vscode/settings.json 中的 cmake.generator 设置

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 源文件
set(SOURCES
    src/atmosphere.cpp
    src/aerodynamics.cpp
    src/propulsion.cpp
    src/trajectory.cpp
    src/integrator.cpp
    src/data_output.cpp
    src/parametric_study.cpp
    src/main.cpp
)

add_executable(cruise_missile_sim ${SOURCES})

target_include_directories(cruise_missile_sim PRIVATE include)

# 输出目录
set_target_properties(cruise_missile_sim PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
