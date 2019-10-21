cmake_minimum_required(VERSION 3.15)
project(HeartValveModel)

set(CMAKE_CXX_STANDARD 17)

add_executable(HeartValveModel main.cpp parameters.h Model.cpp Model.h bezier.cpp bezier.h matplotlibcpp.h)

find_package(PythonLibs 2.7)
target_include_directories(HeartValveModel PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(HeartValveModel ${PYTHON_LIBRARIES})