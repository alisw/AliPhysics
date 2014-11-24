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
include_directories(${AliRoot_SOURCE_DIR}/PHOS/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSbase
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSshuttle
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
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
add_executable(PHOSBCMda PHOSBCMda.cxx) # 1
add_executable(PHOSGAINda PHOSGAINda.cxx) # 2
add_executable(PHOSLEDda PHOSLEDda.cxx) # 3
add_executable(PHOSPEDda PHOSPEDda.cxx) # 4

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(PHOSBCMda PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static ${DATE_MONLIBRARIES} daqDA Root RootExtra ${DATE_RCPROXYLIBRARIES}) # 1
target_link_libraries(PHOSGAINda PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static AmoreDA daqDA ${DATE_MONLIBRARIES}  Root RootExtra ${DATE_RCPROXYLIBRARIES}) # 2
target_link_libraries(PHOSLEDda PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static AmoreDA daqDA ${DATE_MONLIBRARIES} Root RootExtra ${DATE_RCPROXYLIBRARIES}) # 3 
target_link_libraries(PHOSPEDda RAWDatabase-static RAWDatarec-static STEERBase-static AmoreDA daqDA ${DATE_MONLIBRARIES} Root RootExtra ${DATE_RCPROXYLIBRARIES}) # 4

# different flags
set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS} ${AMORE_AUXLIBS}")

set_target_properties(PHOSBCMda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(PHOSBCMda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(PHOSGAINda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(PHOSGAINda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(PHOSLEDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(PHOSLEDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(PHOSPEDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(PHOSPEDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS PHOSBCMda PHOSGAINda PHOSLEDda PHOSPEDda RUNTIME DESTINATION bin)