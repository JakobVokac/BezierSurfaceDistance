# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/s3179222/Desktop/clion_old/clion-2019.2.5/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/s3179222/Desktop/clion_old/clion-2019.2.5/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/s3179222/CLionProjects/HeartValveModel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Parallel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Parallel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Parallel.dir/flags.make

CMakeFiles/Parallel.dir/main.cu.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/main.cu.o: ../main.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/Parallel.dir/main.cu.o"
	/usr/local/cuda/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/s3179222/CLionProjects/HeartValveModel/main.cu -o CMakeFiles/Parallel.dir/main.cu.o

CMakeFiles/Parallel.dir/main.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/Parallel.dir/main.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/Parallel.dir/main.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/Parallel.dir/main.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/Parallel.dir/model.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/model.cpp.o: ../model.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Parallel.dir/model.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/model.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/model.cpp

CMakeFiles/Parallel.dir/model.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/model.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/model.cpp > CMakeFiles/Parallel.dir/model.cpp.i

CMakeFiles/Parallel.dir/model.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/model.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/model.cpp -o CMakeFiles/Parallel.dir/model.cpp.s

CMakeFiles/Parallel.dir/old_code/bezier.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/old_code/bezier.cpp.o: ../old_code/bezier.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Parallel.dir/old_code/bezier.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/old_code/bezier.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/old_code/bezier.cpp

CMakeFiles/Parallel.dir/old_code/bezier.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/old_code/bezier.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/old_code/bezier.cpp > CMakeFiles/Parallel.dir/old_code/bezier.cpp.i

CMakeFiles/Parallel.dir/old_code/bezier.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/old_code/bezier.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/old_code/bezier.cpp -o CMakeFiles/Parallel.dir/old_code/bezier.cpp.s

CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.o: ../old_code/vectorMath.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/old_code/vectorMath.cpp

CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/old_code/vectorMath.cpp > CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.i

CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/old_code/vectorMath.cpp -o CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.s

CMakeFiles/Parallel.dir/old_code/plotting.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/old_code/plotting.cpp.o: ../old_code/plotting.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/Parallel.dir/old_code/plotting.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/old_code/plotting.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/old_code/plotting.cpp

CMakeFiles/Parallel.dir/old_code/plotting.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/old_code/plotting.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/old_code/plotting.cpp > CMakeFiles/Parallel.dir/old_code/plotting.cpp.i

CMakeFiles/Parallel.dir/old_code/plotting.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/old_code/plotting.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/old_code/plotting.cpp -o CMakeFiles/Parallel.dir/old_code/plotting.cpp.s

CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.o: ../geometry/vector/vec3d.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CUDA object CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.o"
	/usr/local/cuda/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/s3179222/CLionProjects/HeartValveModel/geometry/vector/vec3d.cu -o CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.o

CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.o: ../geometry/curve/cubiccrv.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CUDA object CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.o"
	/usr/local/cuda/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/s3179222/CLionProjects/HeartValveModel/geometry/curve/cubiccrv.cu -o CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.o

CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.o: ../geometry/surface/bicubicsrf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/bicubicsrf.cpp

CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/bicubicsrf.cpp > CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.i

CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/bicubicsrf.cpp -o CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.s

CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.o: ../geometry/surface/TopParametric.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/TopParametric.cpp

CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/TopParametric.cpp > CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.i

CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/TopParametric.cpp -o CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.s

CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.o: ../geometry/surface/BottomParametric.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/BottomParametric.cpp

CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/BottomParametric.cpp > CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.i

CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/BottomParametric.cpp -o CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.s

CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.o: ../geometry/surface/mySurface.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/mySurface.cpp

CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/mySurface.cpp > CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.i

CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/mySurface.cpp -o CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.s

CMakeFiles/Parallel.dir/geometry/curve/curve.cu.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/curve/curve.cu.o: ../geometry/curve/curve.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CUDA object CMakeFiles/Parallel.dir/geometry/curve/curve.cu.o"
	/usr/local/cuda/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/s3179222/CLionProjects/HeartValveModel/geometry/curve/curve.cu -o CMakeFiles/Parallel.dir/geometry/curve/curve.cu.o

CMakeFiles/Parallel.dir/geometry/curve/curve.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/Parallel.dir/geometry/curve/curve.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/Parallel.dir/geometry/curve/curve.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/Parallel.dir/geometry/curve/curve.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.o: ../optimizer/optimizer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/optimizer/optimizer.cpp

CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/optimizer/optimizer.cpp > CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.i

CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/optimizer/optimizer.cpp -o CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.s

CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.o: ../optimizer/step/Newton.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/optimizer/step/Newton.cpp

CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/optimizer/step/Newton.cpp > CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.i

CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/optimizer/step/Newton.cpp -o CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.s

CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.o: ../optimizer/preprocessor/bisection.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/optimizer/preprocessor/bisection.cpp

CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/optimizer/preprocessor/bisection.cpp > CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.i

CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/optimizer/preprocessor/bisection.cpp -o CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.s

CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.o: ../optimizer/step/Geometric.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/optimizer/step/Geometric.cpp

CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/optimizer/step/Geometric.cpp > CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.i

CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/optimizer/step/Geometric.cpp -o CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.s

CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.o: ../optimizer/preprocessor/quadraticInterpolation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/optimizer/preprocessor/quadraticInterpolation.cpp

CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/optimizer/preprocessor/quadraticInterpolation.cpp > CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.i

CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/optimizer/preprocessor/quadraticInterpolation.cpp -o CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.s

CMakeFiles/Parallel.dir/measurements/measurements.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/measurements/measurements.cpp.o: ../measurements/measurements.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Building CXX object CMakeFiles/Parallel.dir/measurements/measurements.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/measurements/measurements.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/measurements/measurements.cpp

CMakeFiles/Parallel.dir/measurements/measurements.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/measurements/measurements.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/measurements/measurements.cpp > CMakeFiles/Parallel.dir/measurements/measurements.cpp.i

CMakeFiles/Parallel.dir/measurements/measurements.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/measurements/measurements.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/measurements/measurements.cpp -o CMakeFiles/Parallel.dir/measurements/measurements.cpp.s

CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.o: ../splittingAlgorithm/splittingAlgorithm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_19) "Building CXX object CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/splittingAlgorithm/splittingAlgorithm.cpp

CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/splittingAlgorithm/splittingAlgorithm.cpp > CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.i

CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/splittingAlgorithm/splittingAlgorithm.cpp -o CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.s

CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.o: CMakeFiles/Parallel.dir/flags.make
CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.o: ../geometry/surface/compositeBicubicsrf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_20) "Building CXX object CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.o -c /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/compositeBicubicsrf.cpp

CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/compositeBicubicsrf.cpp > CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.i

CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/s3179222/CLionProjects/HeartValveModel/geometry/surface/compositeBicubicsrf.cpp -o CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.s

# Object files for target Parallel
Parallel_OBJECTS = \
"CMakeFiles/Parallel.dir/main.cu.o" \
"CMakeFiles/Parallel.dir/model.cpp.o" \
"CMakeFiles/Parallel.dir/old_code/bezier.cpp.o" \
"CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.o" \
"CMakeFiles/Parallel.dir/old_code/plotting.cpp.o" \
"CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.o" \
"CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.o" \
"CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.o" \
"CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.o" \
"CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.o" \
"CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.o" \
"CMakeFiles/Parallel.dir/geometry/curve/curve.cu.o" \
"CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.o" \
"CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.o" \
"CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.o" \
"CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.o" \
"CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.o" \
"CMakeFiles/Parallel.dir/measurements/measurements.cpp.o" \
"CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.o" \
"CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.o"

# External object files for target Parallel
Parallel_EXTERNAL_OBJECTS =

Parallel: CMakeFiles/Parallel.dir/main.cu.o
Parallel: CMakeFiles/Parallel.dir/model.cpp.o
Parallel: CMakeFiles/Parallel.dir/old_code/bezier.cpp.o
Parallel: CMakeFiles/Parallel.dir/old_code/vectorMath.cpp.o
Parallel: CMakeFiles/Parallel.dir/old_code/plotting.cpp.o
Parallel: CMakeFiles/Parallel.dir/geometry/vector/vec3d.cu.o
Parallel: CMakeFiles/Parallel.dir/geometry/curve/cubiccrv.cu.o
Parallel: CMakeFiles/Parallel.dir/geometry/surface/bicubicsrf.cpp.o
Parallel: CMakeFiles/Parallel.dir/geometry/surface/TopParametric.cpp.o
Parallel: CMakeFiles/Parallel.dir/geometry/surface/BottomParametric.cpp.o
Parallel: CMakeFiles/Parallel.dir/geometry/surface/mySurface.cpp.o
Parallel: CMakeFiles/Parallel.dir/geometry/curve/curve.cu.o
Parallel: CMakeFiles/Parallel.dir/optimizer/optimizer.cpp.o
Parallel: CMakeFiles/Parallel.dir/optimizer/step/Newton.cpp.o
Parallel: CMakeFiles/Parallel.dir/optimizer/preprocessor/bisection.cpp.o
Parallel: CMakeFiles/Parallel.dir/optimizer/step/Geometric.cpp.o
Parallel: CMakeFiles/Parallel.dir/optimizer/preprocessor/quadraticInterpolation.cpp.o
Parallel: CMakeFiles/Parallel.dir/measurements/measurements.cpp.o
Parallel: CMakeFiles/Parallel.dir/splittingAlgorithm/splittingAlgorithm.cpp.o
Parallel: CMakeFiles/Parallel.dir/geometry/surface/compositeBicubicsrf.cpp.o
Parallel: CMakeFiles/Parallel.dir/build.make
Parallel: /usr/lib/x86_64-linux-gnu/libpython2.7.so
Parallel: CMakeFiles/Parallel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_21) "Linking CXX executable Parallel"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Parallel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Parallel.dir/build: Parallel

.PHONY : CMakeFiles/Parallel.dir/build

CMakeFiles/Parallel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Parallel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Parallel.dir/clean

CMakeFiles/Parallel.dir/depend:
	cd /home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/s3179222/CLionProjects/HeartValveModel /home/s3179222/CLionProjects/HeartValveModel /home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug /home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug /home/s3179222/CLionProjects/HeartValveModel/cmake-build-debug/CMakeFiles/Parallel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Parallel.dir/depend
