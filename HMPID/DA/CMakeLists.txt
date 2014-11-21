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
include_directories(${AliRoot_SOURCE_DIR}/ACORDE/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/HMPID/HMPIDbase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarecOnline
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
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
add_executable(HMPIDda HMPIDda.cxx)
add_executable(HMPIDphysda HMPIDphysda.cxx)

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(HMPIDda HMPIDbase-static RAWDatarecOnline-static RAWDatarec-static RAWDatasim-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(HMPIDphysda HMPIDbase-static RAWDatarecOnline-static RAWDatarec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)

# different flags
set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(HMPIDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(HMPIDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(HMPIDphysda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(HMPIDphysda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS HMPIDda HMPIDphysda RUNTIME DESTINATION bin)