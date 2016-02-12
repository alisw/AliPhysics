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
set(MODULE PWGLFspectra)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/SPECTRA)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGLF/SPECTRA/XtAnalysis
            		    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/JCORRAN/Base
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
		    ${AliPhysics_SOURCE_DIR}/PWG/TRD
  )

# Sources - alphabetical order
set(SRCS
  ChargedHadrons/dNdPt/AlidNdPtUnifiedAnalysisTask.cxx
  ChargedHadrons/dNdPt/MeanPtRaw/AliAnalysisTaskMeanPtRaw.cxx
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
  PiKaPr/ITSsa/AliAnalysisTaskSEITSsaSpectraMultiplicity.cxx
  PiKaPr/TPCTOF/AliAnalysisCombinedHadronSpectra.cxx
  PiKaPr/TPCTOFpA/AliAnalysisTPCTOFpA.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDEvent.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDParticle.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDTrack.cxx
  PiKaPr/TPCTOFfits/AliAnalysisTaskTPCTOFPID.cxx
  PiKaPr/TOF/pp7/TOFSpectrappAnalysis.cxx
  PiKaPr/TOF/pp7/AliAnalysisCombinedHadronSpectra2MC.cxx
  ParticleEfficiency/AliAnalysisTaskParticleEfficiencyLF.cxx
  PiKaPr/TestAOD/AliSpectraAODTrackCuts.cxx
  PiKaPr/TestAOD/AliSpectraAODEventCuts.cxx
  PiKaPr/TestAOD/AliSpectraAODHistoManager.cxx
  PiKaPr/TestAOD/AliSpectraAODPID.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraAOD.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraAllChAOD.cxx
  PiKaPr/TestAOD/AliAnalysisTaskV2AllChAOD.cxx
  Kinks/AliAnalysisKinkESDat.cxx
  Kinks/AliAnalysisKinkESDMC.cxx
  Kinks/AliESDkinkCuts.cxx
  Kinks/AliAnalysisTaskKinksFilimon.cxx
  PiKaPr/TestAOD/AliSpectraBothTrackCuts.cxx
  PiKaPr/TestAOD/AliSpectraBothEventCuts.cxx
  PiKaPr/TestAOD/AliSpectraBothHistoManager.cxx
  PiKaPr/TestAOD/AliSpectraBothPID.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraBoth.cxx
  PiKaPr/TestAOD/AliMultiplictyLoaderTask.cxx
  XtAnalysis/AliXtAnalysis.cxx
  XtAnalysis/AliJXtHistos.cxx
  IdentifiedHighPt/train/AliAnalysisTaskHighPtDeDx.cxx
  IdentifiedHighPt/train/DebugClassesMultESA2013.cxx
  MultEvShape/AliMESeventInfo.cxx
  MultEvShape/AliMEStrackInfo.cxx
  MultEvShape/AliMESbaseTask.cxx
  MultEvShape/AliMEStender.cxx
  MultEvShape/AliMESpidTask.cxx
  MultEvShape/AliMESchgTask.cxx
  MultEvShape/AliMESppColTask.cxx
  AntiprotonToProton/AliProtonAnalysis.cxx
  AntiprotonToProton/AliProtonAnalysisBase.cxx
  AntiprotonToProton/train/AliAnalysisTaskProton.cxx
  AntiprotonToProton/train/AliAnalysisTaskProtonQA.cxx
  AntiprotonToProton/train/AliAnalysisTaskProtonAbsorbtion.cxx
  Spherocity/AliTransverseEventShape.cxx
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Extra files to include in the PARfile
set(PAR_EXTRAFILES
  PiKaPr/TestAOD/HistogramNamesBoth.cxx
  PiKaPr/TestAOD/HistogramNamesBoth.h
  PiKaPr/TestAOD/HistogramNames.cxx
  PiKaPr/TestAOD/HistogramNames.h
  PiKaPr/TestAOD/HistogramsBoth.h
  PiKaPr/TestAOD/Histograms.h
)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependencies
set(LIBDEPS ANALYSISalice PWGTools PWGUDbase EventMixing PWGCFCorrelationsJCORRAN)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "PiKaPr/TestAOD PiKaPr/TPCTOFfits XtAnalysis" "${PAR_EXTRAFILES}")

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
install(DIRECTORY ChargedHadrons/dNdPt/MeanPtRaw/macros DESTINATION
  PWGLF/SPECTRA/ChargedHadrons/dNdPt/MeanPtRaw)
install(DIRECTORY IdentifiedHighPt/train/macros DESTINATION  PWGLF/SPECTRA/IdentifiedHighPt/train)
install(DIRECTORY Kinks/macros DESTINATION  PWGLF/SPECTRA/Kinks)
install(DIRECTORY PiKaPr/TestAOD DESTINATION PWGLF/SPECTRA/PiKaPr)
install(DIRECTORY PiKaPr/ITSsa DESTINATION PWGLF/SPECTRA/PiKaPr)
install(DIRECTORY PiKaPr/TPCTOFfits DESTINATION PWGLF/SPECTRA/PiKaPr)
install(FILES PiKaPr/TPCTOF/AddTaskCombinedHadronSpectra.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOF)
install(FILES PiKaPr/TPCTOFpA/AddTaskTPCTOFpA.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOFpA)
install(DIRECTORY XtAnalysis/macros DESTINATION PWGLF/SPECTRA/XtAnalysis)
install(DIRECTORY MultEvShape DESTINATION PWGLF/SPECTRA FILES_MATCHING PATTERN "*.C")
install(FILES MultEvShape/priorsDist_data_LHC10d_newAliroot.root DESTINATION  PWGLF/SPECTRA/MultEvShape)
install(FILES Spherocity/macros/AddTaskEventShape.C DESTINATION  PWGLF/SPECTRA/Spherocity/macros)
install(DIRECTORY AntiprotonToProton/train DESTINATION PWGLF/SPECTRA/AntiprotonToProton)
