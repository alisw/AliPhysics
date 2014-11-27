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
include_directories(${AliRoot_SOURCE_DIR}/TOF/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/TOF/TOFbase
                   )


# DIM
link_directories(${DIMDIR}/${ODIR})

#daqDA flags
include_directories(${daqDA})
link_directories(${daqDA})

# AMORE definitions
add_definitions(${AMORE_DEFINITIONS})
include_directories(${AMORE_INCLUDE_DIR})

# Generating the DA executable
add_executable(TOFcalibda TOFcalibda.cxx) # 1
add_executable(TOFda TOFda.cxx) # 2
add_executable(TOFnoiseda TOFnoiseda.cxx) # 3
add_executable(TOFphysicsda TOFphysicsda.cxx) # 4
add_executable(TOFpulserda TOFpulserda.cxx) # 5

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(TOFcalibda RAWDatarec-static STEERBase-static TOFbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 1
target_link_libraries(TOFda RAWDatarec-static STEERBase-static TOFbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 2
target_link_libraries(TOFnoiseda RAWDatarec-static STEERBase-static TOFbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 3
target_link_libraries(TOFphysicsda RAWDatarec-static STEERBase-static TOFbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 4
target_link_libraries(TOFpulserda RAWDatarec-static STEERBase-static TOFbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 5

# different flags
set(MODULE_COMPILE_FLAGS "  ${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(TOFcalibda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(TOFcalibda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(TOFda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(TOFda PMDPEDESTALda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(TOFnoiseda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(TOFnoiseda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(TOFphysicsda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(TOFphysicsda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(TOFpulserda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(TOFpulserda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS TOFcalibda
                TOFda
                TOFnoiseda
                TOFphysicsda
                TOFpulserda
        RUNTIME DESTINATION bin)