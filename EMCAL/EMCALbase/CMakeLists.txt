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
set(MODULE EMCALbase)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/EMCAL/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALraw
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliCaloCalibPedestal.cxx
    AliCaloCalibSignal.cxx
    AliEMCALBiasAPD.cxx
    AliEMCALCalibAbs.cxx
    AliEMCALCalibData.cxx
    AliEMCALCalibMapAPD.cxx
    AliEMCALCalibReference.cxx
    AliEMCALCalibTempCoeff.cxx
    AliEMCALCalibTimeDep.cxx
    AliEMCAL.cxx
    AliEMCALDigit.cxx
    AliEMCALDigitizer.cxx
    AliEMCALHit.cxx
    AliEMCALLoader.cxx
    AliEMCALPreprocessor.cxx
    AliEMCALQAChecker.cxx
    AliEMCALRawDigit.cxx
    AliEMCALRawUtils.cxx
    AliEMCALRecParam.cxx
    AliEMCALRecPoint.cxx
    AliEMCALSDigitizer.cxx
    AliEMCALSensorTempArray.cxx
    AliEMCALSensorTemp.cxx
    AliEMCALSimParam.cxx
    AliEMCALSpaceFrame.cxx
    AliEMCALSurvey.cxx
    AliEMCALTrigger.cxx
    AliEMCALTriggerData.cxx
    AliEMCALTriggerDCSConfig.cxx
    AliEMCALTriggerDCSConfigDB.cxx
    AliEMCALTriggerRawDigit.cxx
    AliEMCALTriggerRawDigitMaker.cxx
    AliEMCALTriggerSTUDCSConfig.cxx
    AliEMCALTriggerSTURawStream.cxx
    AliEMCALTriggerTRUDCSConfig.cxx
    AliEMCALWsuCosmicRaySetUp.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER CDB ESD RAWDatabase RAWDatarec EMCALraw EMCALUtils)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(EMCALbase-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library(EMCALbase SHARED $<TARGET_OBJECTS:EMCALbase-object>)

# dependencies
target_link_libraries(${MODULE} STEERBase STEER CDB ESD RAWDatabase RAWDatarec RAWDatasim EMCALraw EMCALUtils)

# Additional compilation flags
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "-O -g -fPIC")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# Static version if DA enabled
if(ALIROOT_STATIC)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:EMCALbase-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})
    target_link_libraries(${MODULE}-static STEERBase-static STEER-static ESD-static CDB-static AOD-static RAWDatabase-static RAWDatarec-static RAWDatasim-static  EMCALraw-static EMCALUtils-static Root)
    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)
