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
include_directories(${AliRoot_SOURCE_DIR}/TRD/DA)

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
add_executable(TRDHALFCHAMBERSTATUSda TRDHALFCHAMBERSTATUSda.cxx) #1
add_executable(TRDPEDESTALda TRDPEDESTALda.cxx) #2
add_executable(TRDVDRIFTda TRDVDRIFTda.cxx) #3

target_link_libraries(TRDHALFCHAMBERSTATUSda CDB-static RAWDatarec-static STEERBase-static TRDbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 1
target_link_libraries(TRDPEDESTALda CDB-static RAWDatarec-static STEERBase-static TRDbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 2
target_link_libraries(TRDVDRIFTda RAWDatarec-static STEERBase-static TRDbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 3

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# different flags
set(MODULE_COMPILE_FLAGS "-O -g ${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(TRDHALFCHAMBERSTATUSda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 1
set_target_properties(TRDHALFCHAMBERSTATUSda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 1

set_target_properties(TRDPEDESTALda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 2
set_target_properties(TRDPEDESTALda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 2

set_target_properties(TRDVDRIFTda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 3
set_target_properties(TRDVDRIFTda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 3

# Installation
install(TARGETS 
                TRDHALFCHAMBERSTATUSda
                TRDPEDESTALda
                TRDVDRIFTda
        RUNTIME DESTINATION bin)