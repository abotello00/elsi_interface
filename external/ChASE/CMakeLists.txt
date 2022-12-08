message("hello from ChASE")

include( FortranCInterface )
FortranCInterface_HEADER( ${CMAKE_BINARY_DIR}/fortran_mangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOL_NAMESPACE "FC_"
  )

#option( CHASE_OUTPUT "ChASE will provide output at each iteration" ON)

list(APPEND includePath ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src)

LIST(APPEND chase_src
  src/interface/chase_c.cpp
  src/interface/chase_f.f90
)

ADD_LIBRARY(chase_f ${chase_src})
set_target_properties(
  "chase_f"
  PROPERTIES Fortran_PREPROCESS ON
  POSITION_INDEPENDENT_CODE TRUE
)

#if( CHASE_OUTPUT )
#  target_compile_definitions( chase_f  INTERFACE "-DCHASE_OUTPUT" )
#endif()

TARGET_INCLUDE_DIRECTORIES(chase_f PRIVATE
  "$<BUILD_INTERFACE:${includePath}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

INSTALL(TARGETS chase_f 
  EXPORT elsiConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )




