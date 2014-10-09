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
set(MODULE dateStream)

include_directories(${CMAKE_SOURCE_DIR}/${MODULE})

# Aditional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/RAW/RAWDatabase
                    ${CMAKE_SOURCE_DIR}/RAW/RAWDatarec
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase
                   )


link_directories(${ROOT_LIBDIR})

# Aditional ROOT libraries that are not shoed by root-config --libs
set(ROOTEXTRALIBS Geom VMC Minuit)

# Add an executable to the project using the specified source files
add_executable(${MODULE} dateStream.cxx)
target_link_libraries(${MODULE} STEERBase RAWDatabase ${ROOTEXTRALIBS} ${ROOT_LIBRARIES})

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)
