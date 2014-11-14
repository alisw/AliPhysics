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

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/alionlinereco
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/MONITOR
                    ${AliRoot_SOURCE_DIR}/MONITOR/MONITORzmq
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

add_definitions(-DZMQ)

if(DIMDIR AND ODIR)
    link_directories(${DIMDIR}/${ODIR})
endif()

if(DATE_FOUND)
    message(STATUS  "Do something here")
endif()

# alistorage executable
add_executable(alionlinereco onlinereco.cxx)
target_link_libraries(alionlinereco MONITOR MONITORzmq STEERBase STEER)
if(PYTHIA6)
    target_link_libraries(alionlinereco EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS alionlinereco RUNTIME DESTINATION bin)


