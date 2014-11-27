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
include_directories(${AliRoot_SOURCE_DIR}/ZDC/DA)

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
add_executable(ZDCCALIBMBda ZDCCALIBMBda.cxx) #1
add_executable(ZDCEMDda ZDCEMDda.cxx) #2
add_executable(ZDCLASERda ZDCLASERda.cxx) #3
add_executable(ZDCMAPPINGda ZDCMAPPINGda.cxx) #4
add_executable(ZDCPEDESTALda ZDCPEDESTALda.cxx) #5

target_link_libraries(ZDCCALIBMBda RAWDatabase-static RAWDatarec-static ZDCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 1
target_link_libraries(ZDCEMDda RAWDatabase-static RAWDatarec-static ZDCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 2
target_link_libraries(ZDCLASERda RAWDatabase-static RAWDatarec-static ZDCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 3
target_link_libraries(ZDCMAPPINGda RAWDatabase-static RAWDatarec-static ZDCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 4
target_link_libraries(ZDCPEDESTALda RAWDatabase-static RAWDatarec-static ZDCbase-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 5

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# different flags
set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(ZDCCALIBMBda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 1
set_target_properties(ZDCCALIBMBda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 1

set_target_properties(ZDCEMDda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 2
set_target_properties(ZDCEMDda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 2
set_target_properties(ZDCLASERda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 3
set_target_properties(ZDCLASERda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 3

set_target_properties(ZDCMAPPINGda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 4
set_target_properties(ZDCMAPPINGda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 4

set_target_properties(ZDCPEDESTALda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 5
set_target_properties(ZDCPEDESTALda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 5

# Installation
install(TARGETS 
                ZDCCALIBMBda
                ZDCEMDda
                ZDCLASERda
                ZDCMAPPINGda
                ZDCPEDESTALda
        RUNTIME DESTINATION bin)