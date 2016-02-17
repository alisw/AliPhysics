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
set(MODULE ITSbase)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/ITS/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/RAW/dateStream
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                   )

# Sources in alphabetical order
set(SRCS
    AliITSBadChannelsSSD.cxx
    AliITSBadChannelsSSDv2.cxx
    AliITSCalibration.cxx
    AliITSCalibrationSDD.cxx
    AliITSCalibrationSPD.cxx
    AliITSCalibrationSSD.cxx
    AliITSChannelDaSSD.cxx
    AliITSChannelStatus.cxx
    AliITSClusterFinder.cxx
    AliITSClusterFinderSDDfast.cxx
    AliITSClusterFinderV2SDD.cxx
    AliITSClusterFinderV2SPD.cxx
    AliITSClusterFinderV2SSD.cxx
    AliITSClusterParam.cxx
    AliITSclusterTable.cxx
    AliITSclusterV2.cxx
    AliITSCompressRawDataSDD.cxx
    AliITSCorrMap1DSDD.cxx
    AliITSCorrMap2DSDD.cxx
    AliITSCorrMapSDD.cxx
    AliITSDCSAnalyzerSDD.cxx
    AliITSDCSDataSDD.cxx
    AliITSDDLModuleMapSDD.cxx
    AliITSDDLRawData.cxx
    AliITSdEdxSamples.cxx
    AliITSDetTypeRec.cxx
    AliITSdigit.cxx
    AliITSdigitSDD.cxx
    AliITSdigitSPD.cxx
    AliITSdigitSSD.cxx
    AliITSDriftSpeedArraySDD.cxx
    AliITSDriftSpeedSDD.cxx
    AliITSEventHeader.cxx
    AliITSFOSignalsSPD.cxx
    AliITSGainSSD.cxx
    AliITSGainSSDv2.cxx
    AliITSgeom.cxx
    AliITSgeomMatrix.cxx
    AliITSgeomTGeo.cxx
    AliITSHLTforSDD.cxx
    AliITSInitGeometry.cxx
    AliITSIntMap.cxx
    AliITSIntMapNode.cxx
    AliITSLoader.cxx
    AliITSMapA1.cxx
    AliITSMapA2.cxx
    AliITSMap.cxx
    AliITSMapSDD.cxx
    AliITSMisAligner.cxx
    AliITSModuleDaSSD.cxx
    AliITSMultReconstructor.cxx
    AliITSNoiseSSD.cxx
    AliITSNoiseSSDv2.cxx
    AliITSPedestalSSD.cxx
    AliITSPedestalSSDv2.cxx
    AliITSPlaneEff.cxx
    AliITSPlaneEffSDD.cxx
    AliITSPlaneEffSPD.cxx
    AliITSPlaneEffSSD.cxx
    AliITSpList.cxx
    AliITSpListItem.cxx
    AliITSQAChecker.cxx
    AliITSQADataMakerRec.cxx
    AliITSQASDDChecker.cxx
    AliITSQASPDChecker.cxx
    AliITSQASSDChecker.cxx
    AliITSQASDDDataMakerRec.cxx
    AliITSQASPDDataMakerRec.cxx
    AliITSQASSDDataMakerRec.cxx
    AliITSRawStream.cxx
    AliITSRawStreamSDDCompressed.cxx
    AliITSRawStreamSDD.cxx
    AliITSRawStreamSPD.cxx
    AliITSRawStreamSPDErrorLog.cxx
    AliITSRawStreamSSD.cxx
    AliITSRawStreamSSDv1.cxx
    AliITSReconstructor.cxx
    AliITSRecoParam.cxx
    AliITSRecPoint.cxx
    AliITSRecPointContainer.cxx
    AliITSresponse.cxx
    AliITSresponseSDD.cxx
    AliITSresponseSPD.cxx
    AliITSsegmentation.cxx
    AliITSsegmentationSDD.cxx
    AliITSsegmentationSPD.cxx
    AliITSsegmentationSSD.cxx
    AliITSSortTrkl.cxx
    AliITSSurveyToAlign.cxx
    AliITStrackMI.cxx
    AliITStrackerMI.cxx
    AliITSTracklPairs.cxx
    AliITStrackSA.cxx
    AliITStrackV2.cxx
    AliITStrackerSA.cxx
    AliITStrackerV2.cxx
    AliITSTransientDigit.cxx
    AliITSTriggerAlgorithmConditions.cxx
    AliITSTriggerConditions.cxx
    AliITSV0Finder.cxx
    AliITSVertexer.cxx
    AliITSVertexer3D.cxx
    AliITSVertexerCosmics.cxx
    AliITSVertexerFast.cxx
    AliITSVertexerFixed.cxx
    AliITSVertexerZ.cxx
    AliITSVertexerZD.cxx
    AliITSZPoint.cxx
   )
if (DATE_FOUND)
    set (SRCS ${SRCS} "AliITSHandleDaSSD.cxx")
    # set the definitions per folder to be passed to ROOT
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS} ${DATE_DYNMONLIBS}")
endif (DATE_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${MODULE_COMPILE_FLAGS}")

# Setting the dependencies
set(ROOT_DEPENDENCIES Core EG GenVector Geom Gpad Graf3d Graf Gui HistPainter Hist MathCore Matrix Minuit Physics RIO Tree)
set(ALIROOT_DEPENDENCIES CDB ESD STEERBase STEER RAWDatabase RAWDatarec RAWDatasim)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)

target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Aditional compile and link flags
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Setting the correct headers for the object as gathered from the dependencies
# After all compilation flags are set
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# Static version if DA enabled
if(ALIROOT_STATIC)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:${MODULE}-object>)
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
