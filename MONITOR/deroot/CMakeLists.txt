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

set(MODULE deroot)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/${MODULE}
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    deroot.cxx
    root2date.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Add an executable to the project using the specified source files
add_executable(${MODULE} ${SRCS})
target_link_libraries(${MODULE} STEERBase RAWDatabase RAWDatarec EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Physics Tree VMC ESD Thread Core RIO Net)

if(ZEROMQ_FOUND)
    target_link_libraries(${MODULE} MONITORzmq)
    add_definitions(-DZMQ)
ENDIF(ZEROMQ_FOUND)

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)
