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

# alieve executable

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/ALIROOT)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_BINARY_DIR}/version
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# aliroot executable
# Sources in alphabetical order
set(SRCS
  aliroot.cxx
)
  
set(CSRCS
  fpu.c
)

add_executable(aliroot ${SRCS} ${CSRCS})
if(ROOT_VERSION_MAJOR LESS 6)
  # ROOT 5: many dependencies (brought by STEER)
  target_link_libraries(aliroot STEER Core Rint)
else()
  # ROOT 6+: minimal dependencies
  target_link_libraries(aliroot STEERBase Core Rint)
endif()
if(TARGET microcern)
  target_link_libraries(aliroot microcern)
endif()

# Installation
install(TARGETS aliroot RUNTIME DESTINATION bin)

# The VDT wrapper
include_directories(SYSTEM ${ROOT_INCLUDE_DIR})
if(EXISTS ${ROOT_INCLUDE_DIR}/vdt/)
  message(STATUS "VDT include directory exists inside ROOT")
  add_library(alivdtwrapper SHARED alivdtwrapper.cxx)
  target_link_libraries(alivdtwrapper dl m)
  install(TARGETS alivdtwrapper DESTINATION lib)
else()
  message(STATUS "VDT not found inside ROOT")
endif()
