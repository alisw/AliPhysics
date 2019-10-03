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

# Module name -> library name
set(MODULE PWGmftmuondepCharmonium)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/mftmuondep/charmonium)

# Aditional include directories in alphabetical order
include_directories(${AliPhysics_SOURCE_DIR}/EVENTMIX
                   )
# Sources - alphabetical order
set(SRCS
    AliAnalysisTaskDimuonBackground.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
# It expect the name of the LinkDef file to be ${MODULE}LinkDef.h
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# VERY IMPORTANT - ADD the correct dependecies
# Seach any Root class you use inside the ROOT rootmaps and add Root library to the dependency
# AliRoot classes are usually located in the name of the include
# Ex: ${AliPhysics_SOURCE_DIR}/ANALYSIS/ANALYSIS -> ANALYSIS
#     ${AliPhysics_SOURCE_DIR}/EVENTMIX - EventMixing
set(ROOT_DEPENDENCIES Core EG Hist MathCore)
set(ALIROOT_DEPENDENCIES STEERBase AOD ANALYSIS ANALYSISalice EventMixing MUONrec MFTbase)

# Generate the ROOT map used to automatically load the dependecies
# Load only the upper level libraries, there is no need to load the dependecies
# If symbols are not found it means that the dependecies are not correctly defined
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Creating the actual library
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${LIBDEPS})

# Custom compilation and linking flags if needed
set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# Setting custom compilation flags for the object
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting customlinking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Library installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Headers installation
install(FILES ${HDRS} DESTINATION include)

# Install the macros
install(FILES AddTaskDimuonBackground.C DESTINATION PWG/mftmuondep/charmonium)

message(STATUS "PWG mftmuondep charmonium enabled")
