cmake_minimum_required(VERSION 3.15)
project(HeartValveModel)

set(CMAKE_CXX_STANDARD 17)

add_executable(HeartValveModel main.cpp parameters.h model.cpp model.h matplotlibcpp.h plotting/plotting.cpp plotting/plotting.h geometry/vector/vec3d.cpp geometry/vector/vec3d.h geometry/curve/cubiccrv.cpp geometry/curve/cubiccrv.h geometry/surface/bicubicsrf.cpp geometry/surface/bicubicsrf.h geometry/surface/surface.h geometry/surface/TopParametric.cpp geometry/surface/TopParametric.h geometry/surface/BottomParametric.cpp geometry/surface/BottomParametric.h geometry/surface/surface.cpp geometry/curve/curve.cpp geometry/curve/curve.h optimizer/optimizer.cpp optimizer/optimizer.h optimizer/step/Newton.cpp optimizer/step/Newton.h optimizer/preprocessor/preprocessor.h optimizer/states.h optimizer/preprocessor/bisection.cpp optimizer/preprocessor/bisection.h optimizer/step/step.h optimizer/step/Geometric.cpp optimizer/step/Geometric.h optimizer/preprocessor/quadraticInterpolation.cpp optimizer/preprocessor/quadraticInterpolation.h measurements/measurements.cpp measurements/measurements.h splittingAlgorithm/splittingAlgorithm.cpp splittingAlgorithm/splittingAlgorithm.h geometry/surface/compositeBicubicsrf.cpp geometry/surface/compositeBicubicsrf.h)

find_package(PythonLibs 2.7)
target_include_directories(HeartValveModel PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(HeartValveModel ${PYTHON_LIBRARIES})