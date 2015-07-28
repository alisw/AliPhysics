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

set(MODULE monitorCheck)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/${MODULE}
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                   )

# Sources in alphabetical order
set(SRCS
    monitorCheck.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Add an executable to the project using the specified source files
add_executable(${MODULE} ${SRCS})
target_link_libraries(${MODULE} RAWDatarec MathCore Thread Graf3d Graf Hist Gpad EG Geom Matrix Physics Minuit VMC Core)

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS} ${DATE_DYNMONLIBS}")
endif(DATE_FOUND)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting target compile and linking flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

if(ZEROMQ_FOUND)
    target_link_libraries(${MODULE} MONITORzmq)
    add_definitions(-DZMQ)
ENDIF(ZEROMQ_FOUND)

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)
