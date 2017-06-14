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
set(MODULE PWGGAPHOSTasks)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/AutoTrendQA
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/CaloCellQA
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/omega3pi
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/ClusterSelection
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/CPVperf
		    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOSCalib
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_Correlation
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_embedding
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_Run2embedding
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_EpRatio
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_GAFlow
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_NeutralMeson
		    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPb
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPb_MC
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPbQA
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_pp_pi0
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_Tagging
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_TriggerQA
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/UserTasks
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Base
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Tasks
                  )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/TENDER/Tender
                    ${AliPhysics_SOURCE_DIR}/TENDER/TenderSupplies
                    ${AliPhysics_SOURCE_DIR}/OADB
		    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PWGGAUtils
                    ${AliPhysics_SOURCE_DIR}/PWGGA/GammaConv/
                    ${AliPhysics_SOURCE_DIR}/PWGCF/FEMTOSCOPY/AliFemto
                   )

# Sources in alphabetical order
set(SRCS
    CaloCellQA/AliAnalysisTaskCaloCellsQA.cxx
    CaloCellQA/AliCaloCellsQA.cxx
    ClusterSelection/AliPHOSTriggerUtils.cxx
    CPVperf/AliAnalysisTaskCPV.cxx
    omega3pi/AliAnalysisTaskOmegaPi0PiPi.cxx
    PHOSCalib/AliAnalysisTaskEmeanCalib.cxx
    PHOSCalib/AliAnalysisTaskPHOSTimeCalib.cxx
    PHOS_Correlation/AliPHOSCorrelations.cxx
    PHOS_embedding/AliAnalysisTaskPi0DiffEfficiency.cxx
    PHOS_embedding/AliAnalysisTaskPi0Efficiency.cxx
    PHOS_embedding/AliPHOSEmbedding.cxx
    PHOS_Run2embedding/AliPHOSEmbeddingRun2.cxx
    PHOS_Run2embedding/AliPHOSEventCounter.cxx
    PHOS_EpRatio/AliAnalysisTaskEpRatio.cxx
    PHOS_GAFlow/AliAnalysisTaskThermalGAFlow.cxx
    PHOS_GAFlow/AliAnalysisTaskThermalGAFlowMC.cxx
    PHOS_NeutralMeson/AliAnalysisTaskPHOSNeutralMeson.cxx
    PHOS_PbPb/AliAnalysisTaskPi0Flow.cxx
    PHOS_PbPb/AliAnalysisTaskGammaFlow.cxx
    PHOS_PbPb/AliAnalysisTaskgg.cxx
    PHOS_PbPb/AliAnalysisTaskggMC.cxx
    PHOS_PbPb/AliAnalysisTaskPHOSPCMgg.cxx
    PHOS_PbPb/AliAnalysisTaskEtaPhigg.cxx
    PHOS_PbPb/AliAnalysisTaskEtaPhiMultgg.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMCAOD.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMC.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMCHijing.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMCParamWeights.cxx
    PHOS_PbPb/AliPHOSTenderTask.cxx
    PHOS_PbPb_MC/AliPHOSHijingEfficiency.cxx
    PHOS_PbPbQA/AliAnalysisTaskPHOSPbPbQA.cxx
    PHOS_pp_pi0/AliAnalysisTaskPi0.cxx
    PHOS_pp_pi0/AliAnalysisTaskPi0Conversion.cxx
    PHOS_pp_pi0/AliAnalysisTaskPi0PP.cxx
    PHOS_pp_pi0/AliCaloPhoton.cxx
    PHOS_Tagging/AliAnalysisTaskTaggedPhotons.cxx
    PHOS_TriggerQA/AliAnalysisTaskPHOSTriggerQA.cxx
    UserTasks/AliAnalysisTaskSEPHOSpPbPi0.cxx
    UserTasks/AliCaloClusterInfo.cxx
    UserTasks/AliPHOSpPbPi0Header.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES Core EG Geom Gpad Hist MathCore Physics RIO Tree)
set(ALIROOT_DEPENDENCIES ANALYSIS ANALYSISalice ESDfilter AOD CDB EMCALUtils ESD PHOSbase PHOSrec PHOSUtils PWGGAUtils STEER STEERBase Tender TenderSupplies)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library. Note the extra includes (relative paths)
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "CaloCellQA omega3pi PHOS_embedding PHOS_PbPbQA PHOS_pp_pi0 PHOS_TriggerQA")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
endif(DATE_FOUND)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting compilation flags for the object
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# Install the macros
install(FILES PHOS_PbPb/AddAODPHOSTender.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSGammaFlow.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSgg.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSggMC.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSPCMgg.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSPi0Flow.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSPi0FlowMC.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSEtaPhi.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(FILES PHOS_PbPb/AddTaskPHOSEtaPhiMult.C DESTINATION PWGGA/PHOSTasks/PHOS_PbPb)
install(DIRECTORY PHOSCalib/macros DESTINATION  PWGGA/PHOSTasks/PHOSCalib)
install(DIRECTORY UserTasks/macros DESTINATION  PWGGA/PHOSTasks/UserTasks)
install(DIRECTORY CaloCellQA/macros DESTINATION PWGGA/PHOSTasks/CaloCellQA)
install(DIRECTORY PHOS_PbPbQA/macros DESTINATION PWGGA/PHOSTasks/PHOS_PbPbQA)
install(DIRECTORY PHOS_TriggerQA/macros DESTINATION PWGGA/PHOSTasks/PHOS_TriggerQA)
install(DIRECTORY PHOS_pPb DESTINATION PWGGA/PHOSTasks)
install(DIRECTORY PHOS_Run2embedding/macros DESTINATION PWGGA//PHOSTasks/PHOS_Run2embedding)
install(DIRECTORY CPVperf/macros DESTINATION PWGGA/PHOSTasks/CPVperf)
install(FILES PHOS_pp_pi0/AddTaskPHOSpp.C DESTINATION PWGGA/PHOSTasks/PHOS_pp_pi0)
install(FILES PHOS_pp_pi0/AddTaskPHOSpi0Conversion.C DESTINATION PWGGA/PHOSTasks/PHOS_pp_pi0)
install(FILES PHOSCalib/AddTaskEmeanCalib.C DESTINATION PWGGA/PHOSTasks/PHOSCalib)
install(FILES PHOSCalib/AddTaskPHOSTimeCalib.C DESTINATION PWGGA/PHOSTasks/PHOSCalib)
install(FILES PHOS_Correlation/AddTaskPi0Correlations.C DESTINATION PWGGA/PHOSTasks/PHOS_Correlation)
install(FILES PHOS_EpRatio/AddTaskPHOSEpRatio.C DESTINATION PWGGA/PHOSTasks/PHOS_EpRatio)
install(FILES PHOS_GAFlow/AddTaskThermalGAFlow.C DESTINATION PWGGA/PHOSTasks/PHOS_GAFlow)
install(FILES PHOS_GAFlow/AddTaskThermalGAFlowMC.C DESTINATION PWGGA/PHOSTasks/PHOS_GAFlow)
install(FILES PHOS_NeutralMeson/AddTaskAliAnalysisTaskPHOSNeutralMeson.C DESTINATION PWGGA/PHOSTasks/PHOS_NeutralMeson)
install(FILES PHOS_Tagging/AddTaskPHOSTagging.C DESTINATION PWGGA/PHOSTasks/PHOS_Tagging)
install(FILES AutoTrendQA/MakeTrendingPHOSQA.C DESTINATION PWGGA/PHOSTasks/AutoTrendQA)
install(FILES AutoTrendQA/DrawTrendingPHOSQA.C DESTINATION PWGGA/PHOSTasks/AutoTrendQA)
