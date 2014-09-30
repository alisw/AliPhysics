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

# Module include folders
include_directories(${CMAKE_SOURCE_DIR}/STEER/${MODULE})

# Additional include directories - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${CMAKE_SOURCE_DIR}/HLT/sim
                    ${CMAKE_SOURCE_DIR}/RAW/
                    ${CMAKE_SOURCE_DIR}/RAW/RAWDatabase
                    ${CMAKE_SOURCE_DIR}/RAW/RAWDatarec
                    ${CMAKE_SOURCE_DIR}/STEER/CDB
                    ${CMAKE_SOURCE_DIR}/STEER/ESD
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase
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
    AliHLTTestInputHandler.cxx
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

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

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
