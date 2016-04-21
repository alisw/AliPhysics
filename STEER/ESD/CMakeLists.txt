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
set(MODULE ESD)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folders
include_directories(${AliRoot_SOURCE_DIR}/STEER/${MODULE})

# Additional include directories - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
		    ${AliRoot_SOURCE_DIR}/RAW/dateStream)

# Sources - alphabetical order
set(SRCS
    AliCascadeVertexer.cxx
    AliESDACORDE.cxx
    AliESDAD.cxx
    AliESDADfriend.cxx
    AliESDCaloCells.cxx
    AliESDCaloCluster.cxx
    AliESDCaloTrigger.cxx
    AliESDcascade.cxx
    AliESDCosmicTrack.cxx
    AliESD.cxx
    AliESDEvent.cxx
    AliESDFIT.cxx
    AliESDFMD.cxx
    AliESDfriend.cxx
    AliESDfriendTrack.cxx
    AliESDHandler.cxx
    AliESDHeader.cxx
    AliESDHLTDecision.cxx
    AliESDHLTtrack.cxx
    AliESDInputHandler.cxx
    AliESDInputHandlerRP.cxx
    AliESDkink.cxx
    AliESDMultITS.cxx
    AliESDMuonCluster.cxx
    AliESDMuonGlobalTrack.cxx
    AliESDMuonPad.cxx
    AliESDMuonTrack.cxx
    AliESDpid.cxx
    AliESDPmdTrack.cxx
    AliESDRun.cxx
    AliESDTOFCluster.cxx
    AliESDTOFHit.cxx
    AliESDTOFMatch.cxx
    AliESDtrack.cxx
    AliESDTrdTrack.cxx
    AliESDTrdTracklet.cxx
    AliESDTrdTrigger.cxx
    AliESDTZERO.cxx
    AliESDTZEROfriend.cxx
    AliESDUtils.cxx
    AliESDv0.cxx
    AliESDV0Params.cxx
    AliESDVertex.cxx
    AliESDVZERO.cxx
    AliESDVZEROfriend.cxx
    AliESDZDC.cxx
    AliExpression.cxx
    AliFMDFloatMap.cxx
    AliFMDMap.cxx
    AliKalmanTrack.cxx
    AliKFParticleBase.cxx
    AliKFParticle.cxx
    AliKFVertex.cxx
    AliMeanVertex.cxx
    AliMultiplicity.cxx
    AliRawDataErrorLog.cxx
    AliSelector.cxx
    AliStrLine.cxx
    AliTrackerBase.cxx
    AliTrackPointArray.cxx
    AliTriggerBCMask.cxx
    AliTriggerClass.cxx
    AliTriggerCluster.cxx
    AliTriggerConfiguration.cxx
    AliTriggerDescriptor.cxx
    AliTriggerInput.cxx
    AliTriggerInteraction.cxx
    AliTriggerIR.cxx
    AliTriggerPFProtection.cxx
    AliTriggerScalersESD.cxx
    AliTriggerScalersRecordESD.cxx
    AliV0vertexer.cxx
    AliVertex.cxx
    AliVertexerTracks.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("ESD" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG Geom Gpad Graf3d MathCore Matrix Physics RIO Tree)
set(ALIROOT_DEPENDENCIES STEERBase)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("ESD" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)

# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:ESD-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# Additional compilation flags for the object
set_target_properties(ESD-object PROPERTIES COMPILE_FLAGS "")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(ESD PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Static version if DA enabled
if(ALIROOT_STATIC)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:ESD-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})

    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")

    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})

    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")
    target_link_libraries(${MODULE}-static ${STATIC_ALIROOT_DEPENDENCIES} Root RootExtra)

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)
