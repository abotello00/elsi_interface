include( FortranCInterface )
FortranCInterface_HEADER( ${CMAKE_BINARY_DIR}/fortran_mangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOL_NAMESPACE "FC_"
  )

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

TARGET_LINK_LIBRARIES(chase_f PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(chase_f PRIVATE
  ${INC_PATHS}	
  "$<BUILD_INTERFACE:${includePath}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

INSTALL(TARGETS chase_f 
  EXPORT elsiConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )




