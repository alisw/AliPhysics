# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************

# Module
set(MODULE microcern)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MICROCERN
                    ${AliRoot_SOURCE_DIR}/MICROCERN/kerngen
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                   )

# Sources in alphabetical order
set(SRCS
   )

# Fortran sources
set(FSRCS
    cltou.F
    gamma.F
    gamma64.F
    sortzv.F
    ufill.F
    fint.F
    abend.F
    cfill.F
    mtlprt.F
    mtlset.F
    kerset.F
    gauss.F
    getenvf.F
    gauss64.F
   )

# C sources
set(CSRCS
    geteni.c
    fchput.c
    fchtak.c
    dummies.c
   )

# On Darwin
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(CSRCS ${CSRCS} lnblnk.c)
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Headers from sources
#string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# System dependent definitions

if(${CMAKE_SYSTEM} MATCHES Darwin)
    add_definitions(-DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC)
endif(${CMAKE_SYSTEM} MATCHES Darwin)
if(${CMAKE_SYSTEM} MATCHES Linux)
    add_definitions(-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ)
endif(${CMAKE_SYSTEM} MATCHES Linux)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
#get_directory_property(incdirs INCLUDE_DIRECTORIES)
#generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")


# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} ${FSRCS} ${CSRCS})
# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

#install(FILES ${HDRS} DESTINATION include)

# Static version
if(ALIROOT_STATIC)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:${MODULE}-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})
    
    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")
    target_link_libraries(${MODULE}-static ${STATIC_ALIROOT_DEPENDENCIES} Root RootExtra)

    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})

    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)
