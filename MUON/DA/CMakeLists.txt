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
include_directories(${AliRoot_SOURCE_DIR}/MUON/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcalib
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/MUON/MUONraw
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
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
add_executable(MUONTRGda MUONTRGda.cxx)
add_executable(MUONTRKGAINda MUONTRKGAINda.cxx)
add_executable(MUONTRKOCCda MUONTRKOCCda.cxx)
add_executable(MUONTRKPEDda MUONTRKPEDda.cxx)

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(MUONTRGda MUONraw MUONcalib MUONraw MUONcore MUONmapping RAWDatarec RAWDatarecOnline ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(MUONTRKGAINda MUONsim MUONshuttle MUONraw MUONmapping MUONcore MUONcalib HLTrec RAWDatarecOnline RAWDatarec ESD STEERBase ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(MUONTRKOCCda MUONraw MUONmapping MUONcore MUONcalib RAWDatarecOnline RAWDatarec RAWDatabase STEERBase ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(MUONTRKPEDda HLTrec MUONsim MUONshuttle MUONraw MUONmapping MUONcore MUONcalib RAWDatarecOnline RAWDatarec ESD STEERBase  ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)

# different flags
set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(MUONTRGda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(MUONTRGda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(MUONTRKGAINda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(MUONTRKGAINda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(MUONTRKOCCda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(MUONTRKOCCda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(MUONTRKPEDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(MUONTRKPEDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS MUONTRGda MUONTRKGAINda MUONTRKOCCda MUONTRKPEDda RUNTIME DESTINATION bin)