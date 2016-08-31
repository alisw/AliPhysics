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
# **************************************************************************/

#Module
set(MODULE EvtGenExternal)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TEvtGen/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/TEvtGen
                    ${AliRoot_SOURCE_DIR}/TEvtGen/EvtGen
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/photosCInterfaces
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/eventRecordInterfaces
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/photosFortranInterfaces
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/utilities)

# Sources - alphabetical order
set(SRCS
  EvtExternalGenFactory.cpp
  EvtExternalGenList.cpp
  EvtPHOTOS.cpp
  EvtPhotosEngine.cpp
  EvtPythia6CommandConverter.cpp
  EvtPythia.cpp
  EvtPythiaEngine.cpp
  EvtTauola.cpp
  EvtTauolaEngine.cpp
 )

# No dictionary

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS})

# Do not generate the ROOT map
# Dependecies
set(LIBDEPS Photos Tauola EvtGen)
if(TARGET pythia8)
  include_directories(${AliRoot_SOURCE_DIR}/PYTHIA8/pythia8/include)
  set(LIBDEPS ${LIBDEPS} pythia8)
  set(EVTGEN_PYTHIA "ON")
  add_definitions("-DEVTGEN_PYTHIA")
endif()

# Additional compilation flags
set(EVTGEN_PHOTOS "ON")
set(EVTGEN_TAUOLA "ON")
add_definitions("-DEVTGEN_PHOTOS -DEVTGEN_TAUOLA")

set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
