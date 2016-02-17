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
set(MODULE Photos)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TEvtGen/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/TEvtGen
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/photosCInterfaces
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/utilities
                    ${AliRoot_SOURCE_DIR}/TEvtGen/Photos/src/photosFortranInterfaces
)

# Sources - alphabetical order
set(SRCS
src/photos-fortran/forW-MEc.cxx
src/photos-fortran/forZ-MEc.cxx
src/photos-fortran/photosC.cxx
src/utilities/Log.cxx
src/utilities/PhotosRandom.cxx
src/utilities/PhotosUtilities.cxx
src/utilities/PhotosDebugRandom.cxx
src/photosFortranInterfaces/PH_HEPEVT_Interface.cxx
src/photosCInterfaces/PhotosBranch.cxx
src/photosCInterfaces/PhotosParticle.cxx
src/photosCInterfaces/Photos.cxx
src/photosCInterfaces/PhotosEvent.cxx
src/eventRecordInterfaces/PhotosHEPEVTEvent.cxx
src/eventRecordInterfaces/PhotosHEPEVTParticle.cxx
src/eventRecordInterfaces/PhotosHepMCEvent.cxx
src/eventRecordInterfaces/PhotosHepMCParticle.cxx
)

# No dictionary

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS})

# Do not generate the ROOT map
# Dependencies
set(LIBDEPS HepMC )

# Additional compilation flags
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
