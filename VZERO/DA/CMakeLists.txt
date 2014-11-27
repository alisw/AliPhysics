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
include_directories(${AliRoot_SOURCE_DIR}/VZERO/DA)

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
add_executable(VZEROCHANNELDELAYTUNINGda VZEROCHANNELDELAYTUNINGda.cxx) #1
add_executable(VZEROda VZEROda.cxx) #2
add_executable(VZEROINTEGRATIONGATETUNINGda VZEROINTEGRATIONGATETUNINGda.cxx) #3
add_executable(VZEROPbPbda VZEROPbPbda.cxx) #4

target_link_libraries(VZEROCHANNELDELAYTUNINGda RAWDatarec-static STEERBase-static VZEROrec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 1
target_link_libraries(VZEROda RAWDatarec-static STEERBase-static VZEROrec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 2
target_link_libraries(VZEROINTEGRATIONGATETUNINGda RAWDatarec-static STEERBase-static VZEROrec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 3
target_link_libraries(VZEROPbPbda RAWDatarec-static STEERBase-static VZEROrec-static ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 4

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# different flags
set(MODULE_COMPILE_FLAGS "  ${DATE_CFLAGS} ${AMORE_CFLAGS}")
set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

set_target_properties(VZEROCHANNELDELAYTUNINGda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 1
set_target_properties(VZEROCHANNELDELAYTUNINGda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 1

set_target_properties(VZEROda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 2
 set_target_properties(VZEROda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 2
 
set_target_properties(VZEROINTEGRATIONGATETUNINGda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 3
set_target_properties(VZEROINTEGRATIONGATETUNINGda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 3

set_target_properties(VZEROPbPbda PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS}) # 4
set_target_properties(VZEROPbPbda PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}") # 4


# Installation
install(TARGETS 
                VZEROCHANNELDELAYTUNINGda
                VZEROda
                VZEROINTEGRATIONGATETUNINGda
                VZEROPbPbda
        RUNTIME DESTINATION bin)