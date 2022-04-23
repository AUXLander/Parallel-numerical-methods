cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

set(PROJECT_NAME LU-decomposition)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} 
                src/main.cpp 
				src/math.h
				src/matrix.h
)
