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
set(MODULE MUONsim)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MUON/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/MUON/MUONbase
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcalib
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON/MUONgeometry
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/MUON/MUONraw
                    ${AliRoot_SOURCE_DIR}/MUON/MUONtrigger
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliMUONChamber.cxx
    AliMUONChamberTrigger.cxx
    AliMUONCommonGeometryBuilder.cxx
    AliMUON.cxx
    AliMUONDigitizerV3.cxx
    AliMUONGainEventGenerator.cxx
    AliMUONHit.cxx
    AliMUONHitStoreV1.cxx
    AliMUONMCDataInterface.cxx
    AliMUONPedestalEventGenerator.cxx
    AliMUONQADataMakerSim.cxx
    AliMUONResponse.cxx
    AliMUONResponseFactory.cxx
    AliMUONResponseTrigger.cxx
    AliMUONResponseTriggerV1.cxx
    AliMUONResponseV0.cxx
    AliMUONSDigitizerV2.cxx
    AliMUONSlatGeometryBuilder.cxx
    AliMUONSt1GeometryBuilder.cxx
    AliMUONSt1GeometryBuilderV2.cxx
    AliMUONSt1SpecialMotif.cxx
    AliMUONSt2GeometryBuilder.cxx
    AliMUONSt2GeometryBuilderV2.cxx
    AliMUONTrigger.cxx
    AliMUONTriggerGeometryBuilder.cxx
    AliMUONv1.cxx
    AliMUONVHitStore.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER CDB ESD RAWDatabase MUONcore MUONbase MUONcalib MUONraw MUONgeometry MUONmapping MUONtrigger)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")


# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set(TARGETFLAGS "-g")

# Optimize
if(CMAKE_COMPILER_IS_GNUCXX)
    if(GCC_MAJOR STREQUAL "4")
        if(GCC_MINOR STREQUAL "3")
            set(TARGETFLAGS "${TARGETFLAGS} -O0")
        else()
            set(TARGETFLAGS "${TARGETFLAGS} -O")
        endif()
    endif()
endif()

set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${TARGETFLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)