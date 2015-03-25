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
include_directories(${AliRoot_SOURCE_DIR}/EVE/alieve_main)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/EVE/EveBase
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# alieve executable
# Sources in alphabetical order
set(SRCS
    alieve_main.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

add_executable(alieve ${SRCS})
target_link_libraries(alieve STEERBase EveBase ASImage)

if(ZEROMQ_FOUND)
    target_link_libraries(alieve MONITORzmq)
    add_definitions(-DZMQ)
endif(ZEROMQ_FOUND)

# alihlteve executable
set(SRCSHLT
    alihlteve_main.cxx
)

add_executable(alihlteve ${SRCSHLT})
target_link_libraries(alihlteve EveBase ASImage)

if(ZEROMQ_FOUND)
    target_link_libraries(alihlteve EveBase MONITORzmq)
endif(ZEROMQ_FOUND)

# Installation
install(TARGETS alieve alihlteve RUNTIME DESTINATION bin)
