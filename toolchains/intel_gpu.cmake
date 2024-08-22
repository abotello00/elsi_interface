### Intel (CUDA) ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fc=ifx -O3 -fp-model precise" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-cc=icx -O3 -fp-model precise -std=c99 -Wno-error=implicit-function-declaration" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-cxx=icpx -O3 -fp-model precise -std=c++11 -Wno-error=implicit-function-declaration" CACHE STRING "C++ flags")
SET(CMAKE_CUDA_FLAGS "-O3 -arch=sm_35" CACHE STRING "CUDA flags")

SET(USE_GPU_CUDA ON CACHE BOOL "Use CUDA-based GPU acceleration in ELPA")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")

SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 $ENV{CUDA_HOME}/lib64" CACHE STRING "External library paths")
SET(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core cusolver cublas cublasLt cudart" CACHE STRING "External libraries")
