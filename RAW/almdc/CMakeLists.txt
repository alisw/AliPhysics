cmake_minimum_required (VERSION 2.8)

# Module
set (MODULE alimdc)
include_directories(${CMAKE_SOURCE_DIR}/${MODULE})
include_directories( ${ROOTSYS}/include ${CMAKE_SOURCE_DIR}/STEER/STEERBase  ${CMAKE_SOURCE_DIR}/RAW/RAWDatabase ${CMAKE_SOURCE_DIR}/RAW/RAWDatarec)

# Standard Root libraries
execute_process(COMMAND root-config --libs OUTPUT_VARIABLE ROOTLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

# Additional Root libraries
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOTLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOTEXTRALIBS "-L${ROOTLIBDIR} -lThread -lGeom -lVMC -lMinuit")

# Add an executable to the project using the specified source files
add_executable(${MODULE} alimdc_main.cxx)
target_link_libraries(${MODULE} STEERBase RAWDatabase ${ROOTEXTRALIBS} ${ROOTLIBS})

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)
