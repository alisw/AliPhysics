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
set(MODULE alirawdump)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/RAW/${MODULE})

# Module aditional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

link_directories(${ROOT_LIBDIR})

# Libraries which are not in the default ROOT list (root-config --libs)
set(ROOTEXTRALIBS Geom VMC Minuit)

# Add an executable to the project using the specified source files
add_executable(${MODULE} alirawdump_main.cxx)

target_link_libraries(${MODULE} STEERBase RAWDatabase EG ESD ${ROOTEXTRALIBS} ${ROOT_LIBRARIES})

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)
