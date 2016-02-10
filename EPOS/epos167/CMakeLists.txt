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
set(MODULE EPOS)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/EPOS
                    ${AliRoot_SOURCE_DIR}/EPOS/epos167
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/EVGEN
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    epostimer.cxx
    ${AliRoot_SOURCE_DIR}/EPOS/TEpos.cxx
    ${AliRoot_SOURCE_DIR}/EPOS/AliGenEpos.cxx
    ${AliRoot_SOURCE_DIR}/EPOS/AliGenEposIsajetToPdgConverter.cxx
   )

# Headers
set(HDRS
    ${AliRoot_SOURCE_DIR}/EPOS/TEpos.h
    ${AliRoot_SOURCE_DIR}/EPOS/AliGenEpos.h
    ${AliRoot_SOURCE_DIR}/EPOS/AliGenEposIsajetToPdgConverter.h
   )

# Fortran sources
set(FSRCS
    epos-app-158.f
    epos-bas-168.f
    epos-con-161.f
    epos-dky-165.f
    epos-dro-168.f
    epos-ems-165.f
    epos-fra-163.f
    epos-ico-158.f
    epos-yyy-168.f
    epos-ids-165.f
    epos-int-168.f
    epos-jps-164.f
    epos-omg-160.f
    epos-par-129.f
    epos-sem-165.f
    epos-rsh-165.f
    epos-qsh-165.f
    epos-tim-155.f
    epos-uti-165.f
    epos-xan-154.f
    epos-xpr-165.f
    eposm.f
    eposu.f
   )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER EVGEN Core EG MathCore)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx ${CSRCS} ${FSRCS})
target_link_libraries(${MODULE} ${LIBDEPS})

# Additional compilation and linking flags
set(CMAKE_C_FLAGS "-ansi ${CMAKE_C_FLAGS}")
set(MODULE_COMPILE_FLAGS "")
set(MODULE_LINK_FLAGS)

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Target properties
set_target_properties(${MODULE} PROPERTIES FORTRAN_FLAGS "")
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Install data files
install(DIRECTORY . DESTINATION EPOS/epos167 FILES_MATCHING PATTERN "*.param")
install(DIRECTORY . DESTINATION EPOS/epos167 FILES_MATCHING PATTERN "*.ini*")
