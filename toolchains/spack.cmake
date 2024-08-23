### GCC ###

SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -ffree-line-length-none -fallow-argument-mismatch" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -std=c++11" CACHE STRING "C++ flags")
SET(OpenMP_Fortran_FLAGS "-fopenmp" CACHE STRING "Fortran OpenMP flags")
SET(OpenMP_Fortran_LIB_NAMES "gomp" CACHE STRING "Fortran OpenMP library name")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_CHASE ON CACHE BOOL "Enable ChASE eigensolver")
SET(CHASE_OUTPUT ON CACHE BOOL "Enable ChASE solver with output")
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_DLAF ON CACHE BOOL "Enable DLA-Future")

SET(LIB_PATHS "/spack-view/lib/" CACHE STRING "External library paths")
SET(INC_PATHS "/spack-view/include/" CACHE STRING "External include paths")
SET(LIBS "libscalapack.so;libopenblas.so;libDLAF_Fortran.so" CACHE STRING "External libraries")

SET(MPIEXEC_1P "mpiexec -n 1 --bind-to core:2" CACHE STRING "Command to run serial tests with 1 MPI task")
SET(MPIEXEC_NP "mpiexec -n 4 --bind-to core:2" CACHE STRING "Command to run parallel tests with multiple MPI tasks")
