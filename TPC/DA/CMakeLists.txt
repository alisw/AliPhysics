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
include_directories(${AliRoot_SOURCE_DIR}/TPC/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
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
add_executable(TPCCEda TPCCEda.cxx) #1
add_executable(TPCCEnewda TPCCEnewda.cxx) #2
add_executable(TPCPEDESTALda TPCPEDESTALda.cxx) #3
add_executable(TPCPULSERda TPCPULSERda.cxx) #4
add_executable(TPCQAda TPCQAda.cxx) #5
add_executable(TPCRAWda TPCRAWda.cxx) #6

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

target_link_libraries(TPCCEda RAWDatarec-static STEERBase-static TPCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} RootExtra Root) # 1
target_link_libraries(TPCCEnewda RAWDatarec-static STEERBase-static TPCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 2
target_link_libraries(TPCPEDESTALda RAWDatarec-static STEERBase-static TPCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 3
target_link_libraries(TPCPULSERda RAWDatarec-static STEERBase-static TPCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 4
target_link_libraries(TPCQAda RAWDatarec-static STEERBase-static TPCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 5
target_link_libraries(TPCRAWda RAWDatarec-static STEERBase-static TPCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 6

# different flags
set(MODULE_COMPILE_FLAGS "-O -g ${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(TPCCEda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 1 
set_target_properties(TPCCEda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 1

set_target_properties(TPCCEnewda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 2 
set_target_properties(TPCCEnewda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 2

set_target_properties(TPCPEDESTALda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 3 
set_target_properties(TPCPEDESTALda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 3

set_target_properties(TPCPULSERda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 4 
set_target_properties(TPCPULSERda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 4

set_target_properties(TPCQAda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 5 
set_target_properties(TPCQAda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 5

set_target_properties(TPCRAWda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 6 
set_target_properties(TPCRAWda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 6


# Installation
install(TARGETS 
                TPCCEda
                TPCCEnewda
                TPCPEDESTALda
                TPCPULSERda
                TPCQAda
                TPCRAWda
        RUNTIME DESTINATION bin)