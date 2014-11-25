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
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                    ${AliRoot_SOURCE_DIR}/ITS/ITSrec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
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
add_executable(ITSSDDBASda ITSSDDBASda.cxx) # 1
add_executable(ITSSDDGAINda ITSSDDGAINda.cxx) # 2
add_executable(ITSSDDINJda ITSSDDINJda.cxx) # 3
add_executable(ITSSPDFOda ITSSPDFOda.cxx) # 4
add_executable(ITSSPDPHYSda ITSSPDPHYSda.cxx) # 5
add_executable(ITSSPDSCANda ITSSPDSCANda.cxx) # 6
add_executable(ITSSPDVertexDiamondda ITSSPDVertexDiamondda.cxx) # 7
add_executable(ITSSSDPEDda ITSSSDPEDda.cxx) # 8

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(ITSSDDBASda ITSbase ITSrec RAWDatarec ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 1
target_link_libraries(ITSSDDGAINda ITSrec-static ITSbase-static RAWDatarec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 2
target_link_libraries(ITSSDDINJda ITSrec-static ITSbase-static RAWDatarec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 3
target_link_libraries(ITSSPDFOda ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 4
target_link_libraries(ITSSPDPHYSda ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 5
target_link_libraries(ITSSPDSCANda ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 6
target_link_libraries(ITSSPDVertexDiamondda ITSrec-static RAWDatarec-static STEERBase-static CDB-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 7
target_link_libraries(ITSSSDPEDda ITSbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 8

# different flags
set(MODULE_COMPILE_FLAGS "-O -g ${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

# 1
set_target_properties(ITSSDDBASda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSDDBASda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 2
set_target_properties(ITSSDDGAINda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSDDGAINda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 3
set_target_properties(ITSSDDINJda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSDDINJda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 4
set_target_properties(ITSSPDFOda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSPDFOda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 5
set_target_properties(ITSSPDPHYSda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSPDPHYSda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 6
set_target_properties(ITSSPDSCANda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSPDSCANda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 7
set_target_properties(ITSSPDVertexDiamondda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSPDVertexDiamondda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# 8
set_target_properties(ITSSSDPEDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
set_target_properties(ITSSSDPEDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")


# Installation
install(TARGETS ITSSDDBASda
                ITSSDDGAINda
                ITSSDDINJda
                ITSSPDFOda
                ITSSPDPHYSda
                ITSSPDSCANda
                ITSSPDVertexDiamondda
                ITSSSDPEDda
                            RUNTIME DESTINATION bin)
