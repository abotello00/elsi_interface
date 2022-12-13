include( FortranCInterface )
FortranCInterface_HEADER( ${CMAKE_BINARY_DIR}/fortran_mangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOL_NAMESPACE "FC_"
  )

list(APPEND includePath ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library( chase INTERFACE )

option( CHASE_OUTPUT "ChASE will provide output at each iteration" OFF)
if( CHASE_OUTPUT )
  add_compile_definitions(CHASE_OUTPUT=1)
endif()

target_include_directories( chase INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>  
)

target_compile_features(chase INTERFACE cxx_auto_type)


IF(USE_GPU_CUDA)
  MESSAGE(STATUS "Enabling GPU (CUDA) support in ChASE")
  enable_language(CUDA)
  # If variable not defined, set it to default value in order to comply with CMake Policy CMP0104
  if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set (CMAKE_CUDA_ARCHITECTURES "80")
  endif()

  # If cuda compute compatibility of your GPU is 8.6 (e.g. RTX 3090) you should build with
  # -DCMAKE_CUDA_ARCHITECTURES=86. In the case you want to build code for more than one 
  # CUDA compute capability (e.g. 70, 75, 80 and 86) then build with
  # -DCMAKE_CUDA_ARCHITECTURES="70;75;80;86"

  # If CMAKE_VERSION < 3.18 then CMake is not compliant with CMAKE policy CMP0104 (introduced 
  # in CMake 3.18) which defines that the variable CMAKE_CUDA_ARCHITECTURES has to be initialized.
  # In that case, the code generation flag has to be set manually.
  # For simplicity and compatibility with newer (3.18+) CMake version, the CMAKE_CUDA_ARCHITECTURES 
  # variable has to be always set, not matter the cmake version.
  if (${CMAKE_VERSION} VERSION_LESS "3.18.0")
    message("CMAKE version less than 3.18. Using old CMAKE CUDA policy for defining CUDA code generation flags")
    foreach(CODE ${CMAKE_CUDA_ARCHITECTURES})
      set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${CODE},code=sm_${CODE} -rdc=true")
    endforeach()
  endif()
  message("-- Configuring the build of CUDA with the compute capability of ${CMAKE_CUDA_ARCHITECTURES}")
  
  add_library( chase_cuda "${CMAKE_CURRENT_LIST_DIR}/src/ChASE-MPI/kernels/shift.cu" )
  set_target_properties( chase_cuda 
	  		 PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

  target_link_libraries( chase_cuda PUBLIC ${LIBS})
ENDIF()

add_library(chase_c src/interface/chase_c.cpp)

IF(USE_GPU_CUDA)
  target_compile_definitions( chase_c PRIVATE HAS_GPU=1)
ENDIF()

TARGET_LINK_LIBRARIES(chase_c INTERFACE chase)
TARGET_LINK_LIBRARIES(chase_c PRIVATE ${LIBS})
IF(USE_GPU_CUDA)
  TARGET_LINK_LIBRARIES(chase_c PRIVATE chase_cuda)
ENDIF()

TARGET_INCLUDE_DIRECTORIES(chase_c PRIVATE
  ${INC_PATHS}
  "$<BUILD_INTERFACE:${includePath}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

ADD_LIBRARY(chase_f src/interface/chase_f.f90)
set_target_properties(
  "chase_f"
  PROPERTIES Fortran_PREPROCESS ON
  POSITION_INDEPENDENT_CODE TRUE
)

TARGET_LINK_LIBRARIES(chase_f INTERFACE chase)
TARGET_LINK_LIBRARIES(chase_f PRIVATE chase_c)
TARGET_INCLUDE_DIRECTORIES(chase_f PRIVATE
  ${INC_PATHS}	
  "$<BUILD_INTERFACE:${includePath}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

IF(USE_GPU_CUDA)
INSTALL(TARGETS chase_f chase chase_c chase_cuda
  EXPORT elsiConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
ELSE()
INSTALL(TARGETS chase_f chase chase_c
  EXPORT elsiConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
ENDIF()
