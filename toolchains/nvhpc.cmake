### NVHPC ###
SET(CMAKE_Fortran_COMPILER mpifort CACHE STRING "MPI Fortran compiler")
set(Fortran_MIN_FLAGS "-O0" CACHE STRING "")
set(FFLAGS "-O0" CACHE STRING "")
SET(CMAKE_C_COMPILER nvc CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER nvc++ CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O2 -Mvect=nosimd" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O1" CACHE STRING "C flags") # SCOTCH segfaults with "O2"
SET(CMAKE_CXX_FLAGS "-O1" CACHE STRING "C++ flags")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_CHASE ON CACHE BOOL "Enable ChASE eigensolver")
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")

SET(LIB_PATHS "$ENV{nvhome}/Linux_x86_64/24.11/comm_libs/12.6/openmpi4/openmpi-4.1.5/lib $ENV{nvhome}/Linux_x86_64/24.11/compilers/lib" CACHE STRING "External library paths")
SET(LIBS "scalapack lapack blas" CACHE STRING "External libraries")

SET(MPIEXEC_1P "mpirun --allow-run-as-root --mca io romio314 -n 1" CACHE STRING "Command to run serial tests with 1 MPI task")
SET(MPIEXEC_NP "mpirun --allow-run-as-root --mca io romio314 -n 4" CACHE STRING "Command to run parallel tests with multiple MPI tasks")
