SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "MPI C++ compiler")
SET(CMAKE_CUDA_HOST_COMPILER "mpicxx" CACHE STRING "CUDA HOST COMPILER")

SET(CMAKE_Fortran_FLAGS "-O3 -ffree-line-length-none -fallow-argument-mismatch" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -std=c++11" CACHE STRING "C++ flags")

SET(CMAKE_CUDA_ARCHITECTURES "60" CACHE STRING "CUDA architecture")
SET(CMAKE_CUDA_FLAGS "-O3 -arch=sm_60" CACHE STRING "CUDA flags")

SET(USE_GPU_CUDA ON CACHE BOOL "Use CUDA-based GPU acceleration in ELPA")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_DLAF ON CACHE BOOL "Enable DLA-Future")

SET(LIB_PATHS "/spack-view/lib/;/spack-view/lib64/" CACHE STRING "External library paths")
SET(INC_PATHS "/spack-view/include/" CACHE STRING "External include paths")
SET(LIBS "libscalapack.so;libopenblas.so;libDLAF_Fortran.so;libcusolver.so;libcublas.so;libcublasLt.so;libcudart.so" CACHE STRING "External libraries")

SET(MPIEXEC_1P "mpiexec -n 1 --bind-to core:2" CACHE STRING "Command to run serial tests with 1 MPI task")
SET(MPIEXEC_NP "mpiexec -n 4 --bind-to core:2" CACHE STRING "Command to run parallel tests with multiple MPI tasks")
