This is the C++ version of the distance computation algorithm for the aortic valve model.

The project includes all of the development, including different types of surfaces and optimization techniques.

The program was made in the CLion IDE so it can be imported as such. It is otherwise configured as a CMake project,
so it can be built or imported as such with other IDEs.

If you have the current version of CMake (3.15) installed on your device, you can simply build an executable of the
current program by running the following in the root directory:

    - make Makefile
    - ./HeartValveModel
    
The program is currently configured to simply take in an input through the **input.txt** file.
The input should be the collection of 3D points you wish to calculate the distance for. 
The input should be just **double** numbers. The algorithm will take each 3 in succession as a 3D point.