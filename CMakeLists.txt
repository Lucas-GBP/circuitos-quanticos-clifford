cmake_minimum_required(VERSION 3.29)
# Nome do projeto
project(CircuitosClifford)
# Define o padrão do compilador C++ (C++11, C++14, etc.)
set(CMAKE_CXX_STANDARD 11)

# Adiciona o executável
add_executable(${PROJECT_NAME}.out src/main.cpp src/libs/qliff.cpp src/utils/circuits.cpp src/utils/bitwise.cpp)