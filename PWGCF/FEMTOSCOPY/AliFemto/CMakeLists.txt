
#Module
set(MODULE PWGCFfemtoscopy)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/FEMTOSCOPY/AliFemto)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
  ${AliPhysics_SOURCE_DIR}/OADB
  ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
  )

# Sources - alphabetical order
set(SRCS
  AliFemtoSimpleAnalysis.cxx
  AliFemtoLikeSignAnalysis.cxx
  AliFemtoVertexAnalysis.cxx
  AliFemtoVertexMultAnalysis.cxx
  AliFemtoAnalysisAzimuthal.cxx
  AliFemtoAnalysisReactionPlane.cxx
  AliFemtoBPLCMS3DCorrFctn.cxx
  AliFemtoBPLCMS3DCorrFctnKK.cxx
  AliFemtoCorrFctn3DLCMSSym.cxx
  AliFemtoQinvCorrFctn.cxx
  AliFemtoBasicEventCut.cxx
  AliFemtoEventCutEstimators.cxx
  AliFemtoBasicTrackCut.cxx
  AliFemtoDummyPairCut.cxx
  AliFemtoCoulomb.cxx
  AliFemtoCutMonitorHandler.cxx
  AliFemtoEvent.cxx
  AliFemtoKink.cxx
  AliFemtoManager.cxx
  AliFemtoPair.cxx
  AliFemtoParticle.cxx
  AliFemtoPicoEvent.cxx
  AliFemtoPicoEventCollectionVectorHideAway.cxx
  AliFemtoTrack.cxx
  AliFemtoV0.cxx
  AliFemtoXi.cxx
  AliFmHelix.cxx
  AliFmPhysicalHelix.cxx
  TpcLocalTransform.cxx
  AliFemtoEventReader.cxx
  AliFemtoEventReaderStandard.cxx
  AliFemtoEventReaderESD.cxx
  AliFemtoEventReaderESDChain.cxx
  AliFemtoEventReaderESDChainKine.cxx
  AliFemtoEventReaderESDChainKineKK.cxx
  AliFemtoModelHiddenInfo.cxx
  AliFemtoModelGlobalHiddenInfo.cxx
  AliFemtoModelGausLCMSFreezeOutGenerator.cxx
  AliFemtoModelWeightGenerator.cxx
  AliFemtoModelWeightGeneratorBasic.cxx
  AliFemtoModelManager.cxx
  AliFemtoModelCorrFctn.cxx
  AliFemtoModelFreezeOutGenerator.cxx
  AliFemtoModelWeightGeneratorLednicky.cxx
  AliFemtoCutMonitorParticleYPt.cxx
  AliFemtoCutMonitorParticleVertPos.cxx
  AliFemtoCutMonitorParticlePID.cxx
  AliFemtoCutMonitorEventMult.cxx
  AliFemtoCutMonitorEventVertex.cxx
  AliFemtoKTPairCut.cxx
  AliFemtoCorrFctnNonIdDR.cxx
  AliFemtoCorrFctn3DSpherical.cxx
  AliFemtoEventReaderAOD.cxx
  AliFemtoEventReaderAODChain.cxx
  AliFemtoEventReaderAODMultSelection.cxx
  AliFemtoAODTrackCut.cxx
  AliFemtoCutMonitor.cxx
  AliFemtoCorrFctn.cxx
  AliAnalysisTaskFemto.cxx
  AliAnalysisTaskFemtoMJ.cxx
  AliTwoTrackRes.cxx
  AliFemtoMCTrackCut.cxx
  AliFemtoEventReaderKinematicsChain.cxx
  AliFemtoEventReaderKinematicsChainESD.cxx
  AliFemtoV0PairCut.cxx
  AliFemtoV0TrackCut.cxx
  AliFemtoCutMonitorV0.cxx
  AliFemtoV0TrackPairCut.cxx
  AliFemtoCutMonitorCollections.cxx
  AliFemtoV0SharedDaughterCut.cxx
  AliFemtoXiSharedDaughterCut.cxx
  AliFemtoSphericityEventCut.cxx
  AliFemtoSpherocityEventCut.cxx
  AliFemtoCorrFctnDPhiStarDEta.cxx
  AliFemtoCorrFctnDPhiStarDEtaStar.cxx
  AliFemtoCorrFctnDPhiStarKStarMergedFraction.cxx
  AliFemtoDeltaPtPairCut.cxx
  AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction.cxx
  AliFemtoXiTrackCut.cxx
  AliFemtoXiPairCut.cxx
  AliFemtoXiTrackPairCut.cxx
  AliFemtoXiV0PairCut.cxx
  AliFemtoBetaTPairCut.cxx
  AliFemtoCutMonitorPairBetaT.cxx
  AliFemtoCutMonitorXi.cxx
  AliFemtoEventReaderAODKinematicsChain.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Additional headers
set(HDRS ${HDRS}
  AliFmThreeVector.h
  AliFmThreeVectorF.h
  AliFmThreeVectorD.h
  AliFmPhysicalHelixD.h
  AliFmHelixD.h
  AliFmLorentzVectorD.h
  AliFemtoV0Cut.h
  AliFemtoKinkCut.h
  AliFemtoXiCut.h
  AliFemtoAnalysisCollection.h
  AliFemtoAnalysis.h
  AliFemtoCorrFctnCollection.h
  AliFemtoEnumeration.h
  AliFemtoHelix.h
  AliFemtoHiddenInfo.h
  AliFemtoLikeSignCorrFctn.h
  AliFemtoLorentzVector.h
  AliFemtoString.h
  AliFemtoThreeVectorD.h
  AliFemtoThreeVector.h
  AliFemtoTypes.h
  AliFemtoVector.h
  AliFmLorentzVector.h
  phys_constants.h
  PhysicalConstants.h
  SystemOfUnits.h
  AliFemtoPairCut.h
  AliFemtoEventCut.h
  AliFemtoParticleCut.h
  AliFemtoTrackCollection.h
  AliFemtoV0Collection.h
  AliFemtoXiCollection.h
  AliFemtoKinkCollection.h
  AliFemtoPicoEventCollection.h
  AliFemtoParticleCollection.h
  AliFemtoCutMonitorCollection.h
  AliFemtoTrackCut.h
  AliFemtoPicoEventCollectionVector.h
  AliFemtoEventWriter.h
  AliFemtoEventWriterCollection.h
  AliFemtoMCTrackCut.h
  AliFemtoEventReaderKinematicsChain.h
  AliFemtoCorrFctnDPhiStarKStarMergedFraction.h
  AliFemtoCorrFctnDPhiStarKStarAverageMergedPointsFraction.h
  AliFemtoBetaTPairCut.h
  AliFemtoCutMonitorPairBetaT.h
  AliFemtoAvgSepCalculator.h
  )

set ( FSRCS
  AliFemtoFsiTools.F
  AliFemtoFsiWeightLednicky.F
  )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx ${FSRCS})

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice OADB)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
