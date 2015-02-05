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

set(MODULE monitorGDC)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/${MODULE}
                   )

# Additional include folders in alphabetical order
include_directories(${ROOT_INCLUDE_DIR}
                              ${AliRoot_SOURCE_DIR}
                   )

if(DATE_FOUND)
    include_directories(
                        ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                        ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                        ${AliRoot_SOURCE_DIR}/STEER/ESD
                        ${AliRoot_SOURCE_DIR}/STEER/STEER
                       )
endif(DATE_FOUND)

# Sources in alphabetical order
set(SRCS
    monitorGDC.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Add an executable to the project using the specified source files
add_executable(${MODULE} ${SRCS})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)
set(MODULE_DEPENDENCIES Core)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS} ${DATE_DYNMONLIBS}")

    set(MODULE_DEPENDENCIES STEER ESD RAWDatarec ITSbase ${MODULE_DEPENDENCIES})
endif(DATE_FOUND)

if(ZEROMQ_FOUND)
    set(MODULE_DEPENDENCIES MONITORzmq ${MODULE_DEPENDENCIES})
    add_definitions(-DZMQ)
ENDIF(ZEROMQ_FOUND)

target_link_libraries(${MODULE} ${MODULE_DEPENDENCIES})

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting target compile and linking flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)
