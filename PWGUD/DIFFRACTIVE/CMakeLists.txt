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
set(MODULE PWGUDdiffractive)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGUD/DIFFRACTIVE
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGUD/DIFFRACTIVE/xsAndTwoProng
                    ${AliPhysics_SOURCE_DIR}/PWGUD/DIFFRACTIVE/legotrain
                    ${AliPhysics_SOURCE_DIR}/PWGUD/DIFFRACTIVE/DGTrigger
                    ${AliPhysics_SOURCE_DIR}/PWGUD/DIFFRACTIVE/CEPdevel
                    ${AliPhysics_SOURCE_DIR}/PWGUD/DIFFRACTIVE/Resonance
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/PWGPP/ITS
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
                   )

# Sources in alphabetical order
set(SRCS
    example/AliAnalysisTaskCDex.cxx
    example/AliCDMesonBaseStripped.cxx
    example/AliCDMesonUtilsStripped.cxx
    legotrain/AliAnalysisTaskCDTree.cxx
    legotrain/AliMultiplicitySelectionCP.cxx
    legotrain/AliAnalysisTaskCDPWA.cxx
    legotrain/AliMultiplicitySelectionCPPWA.cxx
    legotrain/AliAnalysisTaskDG.cxx
    legotrain/AliAnalysisTaskDiffCrossSections.cxx
    legotrain/AliSPDGeometryUtils.cxx
    QA/AliAnalysisTaskCPQA.cxx
    QA/AliAnalysisTaskCEPQA.cxx
    xsAndTwoProng/AliAnalysisTaskCDMeson.cxx
    xsAndTwoProng/AliAnalysisTaskCDskimESD.cxx
    xsAndTwoProng/AliCDMesonBase.cxx
    xsAndTwoProng/AliCDMesonTracks.cxx
    xsAndTwoProng/AliCDMesonUtils.cxx
    DGTrigger/AliTriggerStudy0STP.cxx
    CEPdevel/CEPEventBuffer.cxx
    CEPdevel/CEPTrackBuffer.cxx
    CEPdevel/AliCEPBase.cxx
    CEPdevel/AliCEPUtils.cxx
    CEPdevel/AliAnalysisTaskCEP.cxx
    Resonance/AliAnalysisTaskf0f2.cxx
    Resonance/AliAnalysisPseudoRapidityDensity.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG Geom Hist MathCore Physics RIO Tree)
set(ALIROOT_DEPENDENCIES STEERBase STEER CDB ANALYSISalice ANALYSIS ITSbase ITSrec VZERObase PWGPP)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
endif(DATE_FOUND)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting compilation flags for the object
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# Installing the macros
install(DIRECTORY . DESTINATION PWGUD/DIFFRACTIVE FILES_MATCHING PATTERN "*.C" PATTERN "*.pl" PATTERN "*.sh")
