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
include_directories(${AliRoot_SOURCE_DIR}/HMPID/DA)

# Additional include folders in alphabetical order
include_directories(
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
add_executable(ITSSDDBASda ITSSDDBASda.cxx)
add_executable(ITSSDDGAINda ITSSDDGAINda.cxx)
add_executable(ITSSDDINJda ITSSDDINJda.cxx)
add_executable(ITSSPDFOda ITSSPDFOda.cxx)
add_executable(ITSSPDPHYSda ITSSPDPHYSda.cxx)
add_executable(ITSSPDSCANda ITSSPDSCANda.cxx)
add_executable(ITSSPDVertexDiamondda ITSSPDVertexDiamondda.cxx)
add_executable(ITSSSDPEDda ITSSSDPEDda.cxx)

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(ITSSDDBASda ITSrec-static ITSbase-static RAWDatarec-static RAWDatarecOnline-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(ITSSDDGAINda ITSrec-static ITSbase-static RAWDatarec-static RAWDatarecOnline-static  ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(ITSSDDINJda ITSrec-static ITSbase-static RAWDatarec-static RAWDatarecOnline-static  ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(ITSSPDFOda ITSrec-static ITSbase-static HLTrec-static RAWDatarec-static RAWDatarecOnline-static STEERBase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(ITSSPDPHYSda ITSrec-static ITSbase-static HLTrec-static RAWDatarec-static RAWDatarecOnline-static STEERBase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)
target_link_libraries(ITSSPDSCANda ITSrec-static ITSbase-static HLTrec-static RAWDatarec-static RAWDatarecOnline-static STEERBase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra
target_link_libraries(ITSSPDVertexDiamondda ITSrec-static ITSbase-static HLTrec-static RAWDatarec-static RAWDatarecOnline-static STEERBase-static CDB-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra
target_link_libraries(ITSSDDBASda ITSbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra)

# different flags
set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(HMPIDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(HMPIDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

set_target_properties(HMPIDphysda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(HMPIDphysda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS HMPIDda HMPIDphysda RUNTIME DESTINATION bin)