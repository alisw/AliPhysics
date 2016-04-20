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
set(MODULE STEER)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folders
include_directories(${AliRoot_SOURCE_DIR}/STEER/${MODULE})

# Additional include directories - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${CMAKE_BINARY_DIR}/version/
                   )

# Sources - alphabetical order
set(SRCS
    AliAlignmentTracks.cxx
    AliAlignObj.cxx
    AliAlignObjMatrix.cxx
    AliAlignObjParams.cxx
    AliBaseLoader.cxx
    AliCentralTrigger.cxx
    AliCluster3D.cxx
    AliCluster.cxx
    AliConfig.cxx
    AliCorrQAChecker.cxx
    AliCorrQADataMakerRec.cxx
    AliCTPInputTimeParams.cxx
    AliCTPRawData.cxx
    AliCTPRawStream.cxx
    AliCTPTimeParams.cxx
    AliDataLoader.cxx
    AliDCSArray.cxx
    AliDebugVolume.cxx
    AliDetector.cxx
    AliDetectorRecoParam.cxx
    AliDigit.cxx
    AliDigitizationInput.cxx
    AliDigitizer.cxx
    AliDigitNew.cxx
    AliESDTagCreator.cxx
    AliEventInfo.cxx
    AliFstream.cxx
    AliGausCorr.cxx
    AliGenerator.cxx
    AliGeometry.cxx
    AliGeomManager.cxx
    AliGlobalQAChecker.cxx
    AliGlobalQADataMaker.cxx
    AliGRPManager.cxx
    AliGRPObject.cxx
    AliGRPPreprocessor.cxx
    AliGRPRecoParam.cxx
    AliHelix.cxx
    AliHit.cxx
    AliHitMap.cxx
    AliHLTSimulation.cxx
    AliHLTVEventInputHandler.cxx
    AliKink.cxx
    AliLego.cxx
    AliLegoGenerator.cxx
    AliLegoGeneratorEta.cxx
    AliLegoGeneratorEtaR.cxx
    AliLegoGeneratorPhiZ.cxx
    AliLegoGeneratorXYZ.cxx
    AliLHCClockPhase.cxx
    AliLHCData.cxx
    AliLHCDipValT.cxx
    AliLHCReader.cxx
    AliLoader.cxx
    AliLTUConfig.cxx
    AliMatrixSparse.cxx
    AliMatrixSq.cxx
    AliMC.cxx
    AliMCGenHandler.cxx
    AliMergeCombi.cxx
    AliMillePede2.cxx
    AliMillepede.cxx
    AliMillePedeRecord.cxx
    AliMinResSolve.cxx
    AliMisAligner.cxx
    AliModule.cxx
    AliObjectLoader.cxx
    AliParamList.cxx
    AliParamSolver.cxx
    AliPlaneEff.cxx
    AliQACheckerBase.cxx
    AliQAChecker.cxx
    AliQADataMaker.cxx
    AliQADataMakerRec.cxx
    AliQADataMakerSim.cxx
    AliQAManager.cxx
    AliQAThresholds.cxx
    AliQAv1.cxx
    AliRecoInputHandler.cxx
    AliReconstruction.cxx
    AliReconstructor.cxx
    AliRecoParam.cxx
    AliRecPoint.cxx
    AliRectMatrix.cxx
    AliRelAlignerKalman.cxx
    AliRieman.cxx
    AliRndm.cxx
    AliRun.cxx
    AliRunInfo.cxx
    AliRunLoader.cxx
    AliSelectorRL.cxx
    AliSignalProcesor.cxx
    AliSimulation.cxx
    AliStream.cxx
    AliSurveyObj.cxx
    AliSurveyPoint.cxx
    AliSurveyToAlignObjs.cxx
    AliSymBDMatrix.cxx
    AliSymMatrix.cxx
    AliTracker.cxx
    AliTrackFitter.cxx
    AliTrackFitterKalman.cxx
    AliTrackFitterRieman.cxx
    AliTrackFitterStraight.cxx
    AliTrackleter.cxx
    AliTrackResidualsChi2.cxx
    AliTrackResiduals.cxx
    AliTrackResidualsFast.cxx
    AliTrackResidualsLinear.cxx
    AliTransportMonitor.cxx
    AliTreeLoader.cxx
    AliTriggerDetector.cxx
    AliTriggerRunScalers.cxx
    AliTriggerUtils.cxx
    AliV0.cxx
    AliVectorSparse.cxx
    AliVertexer.cxx
    AliVertexGenerator.cxx
    AliVertexGenFile.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS}
   AliConst.h
   AliRawDataHeaderSim.h
   )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG GenVector GeomPainter Geom Gpad Graf Hist MathCore Matrix Minuit Net Physics Proof RIO Tree VMC)
set(ALIROOT_DEPENDENCIES STEERBase RAWDatabase RAWDatarec CDB ESD ANALYSIS HLTbase)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

set(LIBHLTSIM_VERSION "1")
add_definitions("-DLIBHLTSIM_VERSION=${LIBHLTSIM_VERSION}")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)

target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# Additional compilation flags
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "")

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
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:${MODULE}-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})
    
    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")
    target_link_libraries(${MODULE}-static ${STATIC_ALIROOT_DEPENDENCIES} Root RootExtra)

    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})

    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)

