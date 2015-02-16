# Module
set(MODULE TestShuttle)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/SHUTTLE/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                   )

# DIM include folder and libraries folder
include_directories(${DIM_INCLUDE_DIR})

# Include ApMon.h needed by ROOT monalisa interface
include_directories(${ApMon_INCLUDE_DIR})

# Sources in alphabetical order
set(SRCS
  AliTestDataDCS.cxx
  AliTestPreprocessor.cxx
  AliTestShuttle.cxx
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES CDB STEER)
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES} ${DIM_LIBRARIES})

generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${LIBDEPS})

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)