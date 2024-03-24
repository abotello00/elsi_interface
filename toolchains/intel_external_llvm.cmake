### Intel (External) ###
SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fc=ifx -O3 -ip -fp-model precise -fopenmp" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-cc=icx -O3 -ip -fp-model precise -std=c99 -Wno-error=implicit-function-declaration -fopenmp" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-cxx=icpx -O3 -ip -fp-model precise -std=c++11 -Wno-error=implicit-function-declaration -fopenmp" CACHE STRING "C++ flags")
SET(CMAKE_EXE_LINKER_FLAGS "-lstdc++" CACHE STRING "Linker flags")

# Internal
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_EIGENEXA ON CACHE BOOL "Enable EigenExa")
SET(ENABLE_BSEPACK ON CACHE BOOL "Enable BSEPACK")
SET(ENABLE_CHASE ON CACHE BOOL "Enable ChASE eigensolver")
SET(CHASE_OUTPUT ON CACHE BOOL "Enable ChASE solver with output")

# External
SET(USE_EXTERNAL_ELPA ON CACHE BOOL "Use external ELPA")
SET(USE_EXTERNAL_OMM ON CACHE BOOL "Use external libOMM")
SET(USE_EXTERNAL_PEXSI ON CACHE BOOL "Use external PEXSI")
SET(USE_EXTERNAL_NTPOLY ON CACHE BOOL "Use external NTPoly")
SET(USE_EXTERNAL_BSEPACK ON CACHE BOOL "Use external BSEPACK")

# SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 /home/wy29/opt/elpa/lib /home/wy29/opt/omm/lib /home/wy29/opt/pexsi/lib /home/wy29/opt/SuperLU_DIST_6.4.0/build/SRC /home/wy29/opt/scotch_6.1.1/lib /home/wy29/opt/NTPoly/lib /home/yy244/codes/EigenExa/EigenExa-2.11/install_yy/lib /home/wy29/opt/BSEPACK/lib" CACHE STRING "External library paths")
# SET(INC_PATHS "/home/wy29/opt/elpa/include /home/wy29/opt/omm/include /home/wy29/opt/pexsi/include /home/wy29/opt/SuperLU_DIST_6.4.0/build/SRC /home/wy29/opt/SuperLU_DIST_6.4.0/SRC /home/wy29/opt/scotch_6.1.1/include /home/wy29/opt/NTPoly/include /home/yy244/codes/EigenExa/EigenExa-2.11/install_yy/include" CACHE STRING "External library include paths")
# SET(LIBS "OMM MatrixSwitch elpa pexsi superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotch scotcherr NTPoly EigenExa bsepack sseig mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")

SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 $ENV{ELPA_ROOT}/lib $ENV{OMM_ROOT}/lib $ENV{PEXSI_ROOT}/lib $ENV{SUPERLU_DIST_ROOT}/SRC $ENV{SCOTCH_ROOT}/lib $ENV{NTPOLY_ROOT}/lib $ENV{EIGENEXA_ROOT}/lib $ENV{BSEPACK_ROOT}/lib" CACHE STRING "External library paths")
SET(INC_PATHS "$ENV{ELPA_ROOT}/include $ENV{ELPA_ROOT}/include/modules $ENV{ELPA_ROOT}/include/elpa-2024.03.001.rc1/modules/ $ENV{OMM_ROOT}/include $ENV{PEXSI_ROOT}/include $ENV{PEXSI_ROOT}/../include $ENV{SUPERLU_DIST_ROOT}/../SRC $ENV{SCOTCH_ROOT}/include $ENV{SCOTCH_ROOT}/../include $ENV{NTPOLY_ROOT}/include $ENV{EIGENEXA_ROOT}/include" CACHE STRING "External library include paths")
SET(LIBS "OMM MatrixSwitch elpa pexsi superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotch scotcherr NTPoly EigenExa bsepack sseig mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
