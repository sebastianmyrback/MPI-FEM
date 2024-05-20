### Dependencies
To run this code, you need to have CMake, a compiler that runs C++17, and MPI installed on your system.

### Configure, build, and run
1. Clone repo: ```git clone https://github.com/sebastianmyrback/MPI-FEM.git```.
2. Create build folder, configure and build project:
   ```
   mkdir build
   cd build
   cmake ..
   make
   ```
3. Run main file using mpi with 4 cores: ```mpirun -np 4 ./main```.

The file where the results of my report are taken from is in [examples/main.cpp](examples/main.cpp). Another file, mainly for playing around and computing errors etc is [examples/test.cpp](examples/test.cpp).
