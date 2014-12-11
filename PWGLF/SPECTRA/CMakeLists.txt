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
# **************************************************************************/

#Module
set (MODULE PWGLFspectra)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PWGLF/SPECTRA)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ESDfilter
                    ${AliRoot_SOURCE_DIR}/CORRFW
                    ${AliRoot_SOURCE_DIR}/PWG/Tools
                    ${AliRoot_SOURCE_DIR}/PWGLF/SPECTRA/XtAnalysis
                    ${AliRoot_SOURCE_DIR}/PWGUD/base
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TOF/TOFbase
                    ${AliRoot_SOURCE_DIR}/TOF/TOFrec
  )

# Sources - alphabetical order
set(SRCS
  ChargedHadrons/dNdPt/AlidNdPtAcceptanceCuts.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysis.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysisPbPb.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysisPbPb2011.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysisPbPbAOD.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysispPb.cxx
  ChargedHadrons/dNdPt/AlidNdPtBackgroundCuts.cxx
  ChargedHadrons/dNdPt/AlidNdPtCorrection.cxx
  ChargedHadrons/dNdPt/AlidNdPtCutAnalysis.cxx
  ChargedHadrons/dNdPt/AlidNdPtCutAnalysisPbPb.cxx
  ChargedHadrons/dNdPt/AlidNdPt.cxx
  ChargedHadrons/dNdPt/AlidNdPtEfficiency.cxx
  ChargedHadrons/dNdPt/AlidNdPtEventCuts.cxx
  ChargedHadrons/dNdPt/AlidNdPtHelper.cxx
  ChargedHadrons/dNdPt/AlidNdPtTask.cxx
  ChargedHadrons/dNdPt/AlidNdPtTrackDumpTask.cxx
  ChargedHadrons/dNdPt/AliPtResolAnalysis.cxx
  ChargedHadrons/dNdPt/AliPtResolAnalysisPbPb.cxx
  PiKaPr/ITSsa/AliAnalysisTaskSEITSsaSpectra.cxx
  PiKaPr/TPCTOF/AliAnalysisCombinedHadronSpectra.cxx
  PiKaPr/TPCTOFpA/AliAnalysisTPCTOFpA.cxx
  PiKaPr/TOF/pp7/TOFSpectrappAnalysis.cxx
  PiKaPr/TOF/pp7/AliAnalysisCombinedHadronSpectra2MC.cxx
  ParticleEfficiency/AliAnalysisTaskParticleEfficiency.cxx
  PiKaPr/TestAOD/AliSpectraAODTrackCuts.cxx
  PiKaPr/TestAOD/AliSpectraAODEventCuts.cxx
  PiKaPr/TestAOD/AliSpectraAODHistoManager.cxx
  PiKaPr/TestAOD/AliSpectraAODPID.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraAOD.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraAllChAOD.cxx
  PiKaPr/TestAOD/AliAnalysisTaskV2AllChAOD.cxx
  Kinks/AliAnalysisKinkESDat.cxx
  Kinks/AliAnalysisKinkESDMC.cxx
  PiKaPr/TestAOD/AliSpectraBothTrackCuts.cxx
  PiKaPr/TestAOD/AliSpectraBothEventCuts.cxx
  PiKaPr/TestAOD/AliSpectraBothHistoManager.cxx
  PiKaPr/TestAOD/AliSpectraBothPID.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraBoth.cxx
  Nuclei/B2/AliLnID.cxx
  Nuclei/B2/AliLnHistoMap.cxx
  Nuclei/B2/AliLnAODtrackCuts.cxx
  Nuclei/B2/AliAnalysisTaskB2.cxx
  Nuclei/B2/AliAnalysisTaskB2AOD.cxx
  Nuclei/masses/AliAnalysisNucleiMass.cxx
  Nuclei/masses/AliAnalysisNucleiInfo.cxx
  Nuclei/deuteronpA/AliAnalysisDeuteronpA.cxx
  XtAnalysis/AliXtAnalysis.cxx
  XtAnalysis/AliJXtHistos.cxx
  XtAnalysis/AliJBaseTrack.cxx
  XtAnalysis/AliJBaseCard.cxx
  XtAnalysis/AliJCard.cxx
  IdentifiedHighPt/train/AliAnalysisTaskHighPtDeDx.cxx
  IdentifiedHighPt/train/DebugClassesMultESA2013.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice PWGTools PWGUDbase)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

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

# install macros
install(DIRECTORY ChargedHadrons/dNdPt/macros DESTINATION  PWGLF/SPECTRA/ChargedHadrons/dNdPt)
install(DIRECTORY IdentifiedHighPt/train/macros DESTINATION  PWGLF/SPECTRA/IdentifiedHighPt/train)
install(DIRECTORY Kinks/macros DESTINATION  PWGLF/SPECTRA/Kinks)
install(DIRECTORY Nuclei/B2/macros DESTINATION PWGLF/SPECTRA/Nuclei/B2)
install(DIRECTORY Nuclei/masses  DESTINATION PWGLF/SPECTRA/Nuclei)
install(DIRECTORY PiKaPr/TestAOD DESTINATION PWGLF/SPECTRA/PiKaPr)
install(FILES PiKaPr/TPCTOF/AddTaskCombinedHadronSpectra.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOF)
install(FILES PiKaPr/TPCTOFpA/AddTaskTPCTOFpA.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOFpA)
install(DIRECTORY XtAnalysis/macros DESTINATION PWGLF/SPECTRA/XtAnalysis)