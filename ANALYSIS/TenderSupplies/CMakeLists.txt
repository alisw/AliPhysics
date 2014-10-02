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

# Module
set(MODULE TenderSupplies)

# Module include folder
include_directories(${CMAKE_SOURCE_DIR}/ANALYSIS/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${CMAKE_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${CMAKE_SOURCE_DIR}/ANALYSIS/Tender
                    ${CMAKE_SOURCE_DIR}/EMCAL
                    ${CMAKE_SOURCE_DIR}/EMCAL/EMCALbase
                    ${CMAKE_SOURCE_DIR}/EMCAL/EMCALrec
                    ${CMAKE_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${CMAKE_SOURCE_DIR}/CORRFW
                    ${CMAKE_SOURCE_DIR}/OADB
                    ${CMAKE_SOURCE_DIR}/PHOS
                    ${CMAKE_SOURCE_DIR}/PHOS/PHOSbase
                    ${CMAKE_SOURCE_DIR}/PHOS/PHOSrec
                    ${CMAKE_SOURCE_DIR}/PHOS/PHOSUtils
                    ${CMAKE_SOURCE_DIR}/STEER/AOD
                    ${CMAKE_SOURCE_DIR}/STEER/CDB
                    ${CMAKE_SOURCE_DIR}/STEER/ESD
                    ${CMAKE_SOURCE_DIR}/STEER/STEER
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase
                    ${CMAKE_SOURCE_DIR}/T0
                    ${CMAKE_SOURCE_DIR}/T0/T0rec
                    ${CMAKE_SOURCE_DIR}/TOF/
                    ${CMAKE_SOURCE_DIR}/TOF/TOFbase
                    ${CMAKE_SOURCE_DIR}/TOF/TOFrec
                    ${CMAKE_SOURCE_DIR}/TPC
                    ${CMAKE_SOURCE_DIR}/TPC/Base
                    ${CMAKE_SOURCE_DIR}/TPC/TPCbase
                    ${CMAKE_SOURCE_DIR}/TRD
                    ${CMAKE_SOURCE_DIR}/TRD/Cal
                    ${CMAKE_SOURCE_DIR}/TRD/TRDbase
                    ${CMAKE_SOURCE_DIR}/VZERO
                    ${CMAKE_SOURCE_DIR}/VZERO/VZERObase
                    ${CMAKE_SOURCE_DIR}/VZERO/VZEROrec
                   )

# Sources
set(SRCS
    AliAnalysisTaskVZEROEqFactorTask.cxx
    AliEMCALTenderSupply.cxx
    AliHMPIDTenderSupply.cxx
    AliPHOSTenderSupply.cxx
    AliPIDTenderSupply.cxx
    AliT0TenderSupply.cxx
    AliTOFTenderSupply.cxx
    AliTPCTenderSupply.cxx
    AliTrackFixTenderSupply.cxx
    AliTRDTenderSupply.cxx
    AliVtxTenderSupply.cxx
    AliVZEROTenderSupply.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSIS ANALYSISalice Tender EMCAL EMCALbase EMCALrec EMCALUtils CORRFW OADB PHOS PHOSbase PHOSrec PHOSUtils AOD CDB ESD STEER STEERBase T0 T0rec TOF TOFbase TOFrec TPC TPCBase TRD TRDCal TRDBase VZERO VZEROBase VZEROrec)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "-O -g")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

