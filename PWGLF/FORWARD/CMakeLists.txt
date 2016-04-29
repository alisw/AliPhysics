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
set(MODULE PWGLFforward2)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/FORWARD)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGLF/FORWARD/analysis2
  )

# Sources under analysis2/ - alphabetical order
set(SRCS_analysis2
  # AOD objects
  analysis2/AliAODCentralMult.cxx
  analysis2/AliAODForwardEP.cxx
  analysis2/AliAODForwardMult.cxx
  analysis2/AliAODMultEventClass.cxx
  # Interface to OADB backed by a TTree 
  analysis2/AliOADBForward.cxx
  # Correction managers 
  analysis2/AliCorrectionManagerBase.cxx 
  analysis2/AliForwardCorrectionManager.cxx
  analysis2/AliCentralCorrectionManager.cxx
  # Central (SPD) code 
  analysis2/AliCentralCorrAcceptance.cxx 
  analysis2/AliCentralCorrSecondaryMap.cxx
  analysis2/AliCentraldNdetaTask.cxx
  analysis2/AliCentralMCCorrectionsTask.cxx
  analysis2/AliCentralMCMultiplicityTask.cxx
  analysis2/AliCentralMultiplicityTask.cxx
  analysis2/AliSPDMCTrackDensity.cxx
  # Aux tasks and code 
  analysis2/AliCopyHeaderTask.cxx
  analysis2/AliBasedNdetaTask.cxx
  analysis2/AliBaseMCTrackDensity.cxx
  analysis2/AliMCTruthdNdetaTask.cxx
  analysis2/AliDisplacedVertexSelection.cxx
  analysis2/AliDisplacedVertexSelectionAD.cxx
  analysis2/AliPoissonCalculator.cxx
  analysis2/AliMCAuxHandler.cxx
  analysis2/AliTestAD.cxx
  # Forward Aux
  analysis2/AliForwardUtil.cxx
  analysis2/AliFMDMultCuts.cxx
  # FMD corrections
  analysis2/AliFMDCorrAcceptance.cxx
  analysis2/AliFMDCorrDoubleHit.cxx
  analysis2/AliFMDCorrector.cxx 
  analysis2/AliFMDCorrELossFit.cxx 
  analysis2/AliFMDCorrMergingEfficiency.cxx
  analysis2/AliFMDCorrSecondaryMap.cxx 
  analysis2/AliFMDCorrVertexBias.cxx
  # FMD algortithmns 
  analysis2/AliFMDDensityCalculator.cxx
  analysis2/AliFMDEnergyFitter.cxx
  analysis2/AliFMDEventInspector.cxx
  analysis2/AliFMDEventPlaneFinder.cxx
  analysis2/AliFMDHistCollector.cxx
  analysis2/AliFMDESDFixer.cxx
  analysis2/AliFMDSharingFilter.cxx
  analysis2/AliMultEventClassifier.cxx
  # FMD MC algorithms
  analysis2/AliFMDMCCorrector.cxx
  analysis2/AliFMDMCDensityCalculator.cxx
  analysis2/AliFMDMCEventInspector.cxx
  analysis2/AliFMDMCSharingFilter.cxx
  analysis2/AliFMDMCTrackDensity.cxx
  # Forward tasks 
  analysis2/AliBaseESDTask.cxx
  analysis2/AliBaseAODTask.cxx
  analysis2/AliBaseMCCorrectionsTask.cxx
  analysis2/AliFMDEnergyFitterTask.cxx 
  analysis2/AliFMDEventPlaneTask.cxx
  analysis2/AliForwarddNdetaTask.cxx
  analysis2/AliForwardFlowTaskQC.cxx
  analysis2/AliForwardMCCorrectionsTask.cxx
  analysis2/AliForwardMCFlowTaskQC.cxx
  analysis2/AliForwardMCMultiplicityTask.cxx 
  analysis2/AliForwardMultiplicityBase.cxx
  analysis2/AliForwardMultiplicityTask.cxx
  analysis2/AliForwardQATask.cxx
  analysis2/AliBaseMultTask.cxx
  analysis2/AliForwardTriggerBiasCorrection.cxx
  analysis2/AliForwardCreateResponseMatrices.cxx
  analysis2/AliForwardMultiplicityDistribution.cxx
  analysis2/AliForwardMultDists.cxx
  # MC investigations 
  analysis2/AliFMDMCTrackELoss.cxx
  analysis2/AliFMDMCTrackInspector.cxx 
  analysis2/AliFMDMCTrackInspectorTask.cxx
  # MC Weights
  analysis2/AliBaseMCWeights.cxx	
  analysis2/AliForwardFlowWeights.cxx
  analysis2/AliSimplePidWeights.cxx
  analysis2/AliPtEtaPidWeights.cxx
  # Other
  analysis2/AliMultEventClassifierTask.cxx
)

# Sources under GEO/ - alphabetical order
set(SRCS_GEO
  # ZDC tasks
  GEO/AliAnalysisTaskZDCPbPb.cxx
  GEO/AliAnalysisTaskZDCTreeMaker.cxx
  GEO/AliAnalysisTaskZDCpAcalib.cxx
)

# Sources under photons/ - alphabetical order
set(SRCS_photons
  # Photons 
  photons/AliAnalysisTaskPMD.cxx
  photons/AliAnalysisTaskPMDSim.cxx
)

# Putting all sources together
set(SRCS ${SRCS_analysis2} ${SRCS_GEO} ${SRCS_photons})

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Additional headers for analysis2/
set(HDRS_analysis2
  analysis2/AliFMDStripIndex.h 
  analysis2/AliFMDEncodedEdx.h 
  analysis2/AliLandauGaus.h 
  analysis2/AliLandauGausFitter.h 
  analysis2/AliFMDCorrNoiseGain.h
)

# Add additional headers to the common list of headers
list(APPEND HDRS ${HDRS_analysis2})

# Additional definitions
if((DEFINED ALIROOT_SVN_REVISION)) 
  add_definitions(-DALIROOT_SVN_REVISION=${ALIROOT_SVN_REVISION})
elseif((DEFINED ALIROOT_REVISION)) 
  add_definitions(-DALIROOT_REVISION=\"${ALIROOT_REVISION}\")
endif((DEFINED ALIROOT_SVN_REVISION)) 
if((DEFINED ALIROOT_SVN_BRANCH)) 
  add_definitions(-DALIROOT_SVN_BRANCH=\"${ALIROOT_SVN_BRANCH}\")
elseif((DEFINED ALIROOT_BRANCH)) 
  add_definitions(-DALIROOT_BRANCH=\"${ALIROOT_BRANCH}\")
endif((DEFINED ALIROOT_SVN_BRANCH))

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate a PARfile target for this library
add_target_parfile(${MODULE}
  "${SRCS}"
  "${HDRS}"
  "${MODULE}LinkDef.h"
  "${LIBDEPS}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice OADB)
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

# Extra executable
if(ROOT_HASALIEN AND AliRoot_HASALIEN)
  add_executable(runTrain
                 trains/trainMain.cxx)
  target_link_libraries(runTrain ANALYSIS ANALYSISalice OADB Core Rint RIO)
  install(TARGETS runTrain
          RUNTIME DESTINATION bin)
endif(ROOT_HASALIEN AND AliRoot_HASALIEN)
# --------------------------------------------------------------------
# Extra targets 
#

add_custom_command( OUTPUT doc/html/index.html 
		    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/doc
                    COMMAND sed -e 's,\\.\\./,${CMAKE_CURRENT_SOURCE_DIR}/,' -e 's,OUTPUT_DIRECTORY *=.*,OUTPUT_DIRECTORY = doc/,' < ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile > ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
		    COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
		    DEPENDS ${SRCS} ${HDRS})

add_custom_target( PWGLFforward-doxy DEPENDS doc/html/index.html )

# --------------------------------------------------------------------
# Extra installation targets 
#
# Install scripts 
install ( DIRECTORY analysis2/scripts 
	  DESTINATION PWGLF/FORWARD/analysis2 
	  USE_SOURCE_PERMISSIONS
	  PATTERN "*~"    EXCLUDE)
	  
# install QA scripts 
install ( DIRECTORY analysis2/qa
	  DESTINATION PWGLF/FORWARD/analysis2 
	  USE_SOURCE_PERMISSIONS
	  PATTERN "*~"    EXCLUDE)
# install correction scripts 
install ( DIRECTORY analysis2/corrs
	  DESTINATION PWGLF/FORWARD/analysis2 
	  PATTERN "*~"    EXCLUDE)
# install trains 
install ( DIRECTORY analysis2/trains 
	  DESTINATION PWGLF/FORWARD/analysis2 
	  PATTERN "*~"   EXCLUDE)
# install trains 
install ( DIRECTORY trains 
	  DESTINATION PWGLF/FORWARD 
	  FILES_MATCHING PATTERN "*.C"
	  PATTERN "*~"        EXCLUDE)
# install dN/deta scripts 
install ( DIRECTORY analysis2/dndeta
	  DESTINATION PWGLF/FORWARD/analysis2 
	  USE_SOURCE_PERMISSIONS
	  PATTERN "*~"    EXCLUDE)
# install GraphSysErr scripts 
install ( DIRECTORY analysis2/gse
	  DESTINATION PWGLF/FORWARD/analysis2 
	  USE_SOURCE_PERMISSIONS
	  PATTERN "*~"    EXCLUDE)
# Install other scripts
install ( FILES analysis2/AddTaskdNdeta.C
  		analysis2/AddTaskCentralMCCorr.C
		analysis2/AddTaskCentralMult.C
  		analysis2/AddTaskCentralTracks.C
                analysis2/AddTaskCopyHeader.C
		analysis2/AddTaskFMDELoss.C
		analysis2/AddTaskFMDEventPlane.C
		analysis2/AddTaskFMDMCTracks.C
		analysis2/AddTaskForwardFlowQC.C
  		analysis2/AddTaskForwardMCCorr.C
                analysis2/AddTaskForwardMult.C
                analysis2/AddTaskForwardMultDists.C
                analysis2/AddTaskForwardQA.C
		analysis2/AddTaskMCParticleFilter.C
		analysis2/AddTaskMultDistributions.C
		analysis2/AddTaskCreateRespMatr.C
		analysis2/AddTaskTriggerCorrection.C 
		analysis2/DrawdNdeta.C
		analysis2/DrawdNdetaSummary.C
		analysis2/DrawMCCorrSummary.C
		analysis2/DrawAODSummary.C
		analysis2/DrawMultDistsSummary.C
		analysis2/DrawUnfoldedSummary.C
		analysis2/ForwardAODConfig.C
		analysis2/CentralAODConfig.C
		analysis2/dNdetaConfig.C
		analysis2/elossFitConfig.C
                analysis2/MakeAOD.C
		analysis2/MakedNdeta.C
                analysis2/MakeELossFits.C
		analysis2/MakeEvaluateTriggers.C
		analysis2/MakeFlow.C
		analysis2/MakeFMDEventPlane.C
                analysis2/MakeForwardQA.C
		analysis2/MakeMCCorr.C		
           DESTINATION PWGLF/FORWARD/analysis2 )

# Install corrections 
install ( DIRECTORY corrections/Empirical
		    corrections/Normalization
          DESTINATION PWGLF/FORWARD/corrections 
	  PATTERN "*~"   EXCLUDE
	  PATTERN ".svn" EXCLUDE)

# Install shell script 
install ( FILES analysis2/Run.sh 
		analysis2/liteAnalysis.sh
		analysis2/gridAnalysis.sh
		analysis2/baseAnalysis.sh
	  DESTINATION PWGLF/FORWARD/analysis2 
	  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                      GROUP_EXECUTE GROUP_READ
		      WORLD_EXECUTE WORLD_READ )

# install sim 
install ( DIRECTORY analysis2/sim 
         DESTINATION PWGLF/FORWARD/analysis2)

install(DIRECTORY macros DESTINATION PWGLF/FORWARD)
# --------------------------------------------------------------------

