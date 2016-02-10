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
set(MODULE PHOSbase)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PHOS/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSUtils
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliPHOS.cxx
    AliPHOSCalibData.cxx
    AliPHOSClusterizer.cxx
    AliPHOSClusterizerv1.cxx
    AliPHOSCpv3GConnection.cxx
    AliPHOSCpvBadChannelsMap.cxx
    AliPHOSCpvCalibData.cxx
    AliPHOSCpvGainCalibDA.cxx
    AliPHOSCpvParam.cxx
    AliPHOSCpvPedProducer.cxx
    AliPHOSCpvRawDigiProducer.cxx
    AliPHOSCpvRawStream.cxx
    AliPHOSCpvRecPoint.cxx
    AliPHOSDigit.cxx
    AliPHOSDigitizer.cxx
    AliPHOSEmcBadChannelsMap.cxx
    AliPHOSEmcCalibData.cxx
    AliPHOSEmcRecPoint.cxx
    AliPHOSFastRecParticle.cxx
    AliPHOSGeometry.cxx
    AliPHOSHit.cxx
    AliPHOSImpact.cxx
    AliPHOSLink.cxx
    AliPHOSLoader.cxx
    AliPHOSModuleMisalignment.cxx
    AliPHOSPID.cxx
    AliPHOSPulseGenerator.cxx
    AliPHOSQAChecker.cxx
    AliPHOSRawDigiProducer.cxx
    AliPHOSRawFitterv0.cxx
    AliPHOSRawFitterv1.cxx
    AliPHOSRawFitterv2.cxx
    AliPHOSRawFitterv3.cxx
    AliPHOSRawFitterv4.cxx
    AliPHOSRecoParam.cxx
    AliPHOSReconstructor.cxx
    AliPHOSRecParticle.cxx
    AliPHOSRecPoint.cxx
    AliPHOSPIDv1.cxx
    AliPHOSSDigitizer.cxx
    AliPHOSSimParam.cxx
    AliPHOSSurvey1.cxx
    AliPHOSSurvey.cxx
    AliPHOSTick.cxx
    AliPHOSTracker.cxx
    AliPHOSTrackSegment.cxx
    AliPHOSTrackSegmentMaker.cxx
    AliPHOSTrackSegmentMakerv1.cxx
    AliPHOSTrigger.cxx
    AliPHOSTriggerParameters.cxx
    AliPHOSTriggerRawDigit.cxx
    AliPHOSTriggerRawDigiProducer.cxx
    AliPHOSTriggerRawReader.cxx
    AliPHOSTriggerSTURawStream.cxx
    AliPHOSTRURawReader.cxx

   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG Geom Gpad Graf Hist MathCore Matrix Minuit Physics RIO Tree VMC)
set(ALIROOT_DEPENDENCIES STEERBase STEER CDB ESD RAWDatabase RAWDatarec RAWDatasim PHOSUtils)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Follow headers dependencies
add_dependencies(${MODULE}-object ${ALIROOT_DEPENDENCIES})
# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:PHOSbase-object>)

# Linking
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Additional compilation flags
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS " ")

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
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:PHOSbase-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})

    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")
    target_link_libraries(${MODULE}-static ${STATIC_ALIROOT_DEPENDENCIES} Root)

    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})

    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)
