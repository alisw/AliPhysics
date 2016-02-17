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
set(MODULE HepMC         )
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TEvtGen/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/TEvtGen
                    ${AliRoot_SOURCE_DIR}/TEvtGen/EvtGen/EvtGenModels
)

# Sources - alphabetical order
set(SRCS
  CompareGenEvent.cc
  filterEvent.cc
  Flow.cc
  GenCrossSection.cc
  GenEvent.cc
  GenEventStreamIO.cc
  GenParticle.cc
  GenRanges.cc
  GenVertex.cc
  HeavyIon.cc
  HEPEVT_Wrapper.cc
  HerwigWrapper.cc
  IO_AsciiParticles.cc
  IO_GenEvent.cc
  IO_HEPEVT.cc
  IO_HERWIG.cc
  PdfInfo.cc
  Polarization.cc
  SearchVector.cc
  SimpleVector.icc
  StreamHelpers.cc
  StreamInfo.cc
  WeightContainer.cc
  Units.cc 
)

set (FSRCS
  dummy.F
)

# No dictionary

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} ${FSRCS})

# Do not generate the ROOT map

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")
target_link_libraries(${MODULE} )

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
