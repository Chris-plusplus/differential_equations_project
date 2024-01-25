cmake_minimum_required(VERSION 3.19)

project("differential_equations_project")

set(CMAKE_CXX_STANDARD 20)

if (MSVC)
    add_compile_options("/Zc:__cplusplus")
    add_compile_options("/Zc:preprocessor")
endif (MSVC)

file(GLOB_RECURSE MY_SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} "main.cpp")
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_sources(${PROJECT_NAME} PUBLIC ${MY_SOURCES})