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

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/CaloCellQA
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/omega3pi
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_Correlation
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_embedding
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_EpRatio
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPb
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPb_MC
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPbQA
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_pp_pi0
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_Tagging
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_TriggerQA
                    ${AliRoot_SOURCE_DIR}/PWGGA/PHOSTasks/UserTasks
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/Tender
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/TenderSupplies
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSbase
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSrec
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSUtils
                    ${AliRoot_SOURCE_DIR}/PWGGA/PWGGAUtils
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    CaloCellQA/AliAnalysisTaskCaloCellsQA.cxx
    CaloCellQA/AliCaloCellsQA.cxx
    omega3pi/AliAnalysisTaskOmegaPi0PiPi.cxx
    PHOS_Correlation/AliPHOSCorrelations.cxx
    PHOS_embedding/AliAnalysisTaskPi0DiffEfficiency.cxx
    PHOS_embedding/AliAnalysisTaskPi0Efficiency.cxx
    PHOS_embedding/AliPHOSEmbedding.cxx
    PHOS_EpRatio/AliAnalysisTaskEpRatio.cxx
    PHOS_PbPb/AliAnalysisTaskPi0Flow.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMCAOD.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMC.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMCHijing.cxx
    PHOS_PbPb/AliAnalysisTaskPi0FlowMCParamWeights.cxx
    PHOS_PbPb/AliPHOSTenderTask.cxx
    PHOS_PbPb_MC/AliPHOSHijingEfficiency.cxx
    PHOS_PbPbQA/AliAnalysisTaskPHOSPbPbQA.cxx
    PHOS_pp_pi0/AliAnalysisTaskPi0.cxx
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
set(ALIROOT_DEPENDENCIES ANALYSIS ANALYSISalice AOD CDB EMCALUtils ESD PHOSbase PHOSrec PHOSUtils PWGGAUtils STEER STEERBase Tender TenderSupplies)

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

install(FILES
                CaloCellQA/macros/AddTaskCaloCellsQA.C
                CaloCellQA/macros/AddTaskPHOSQA.C
       DESTINATION
                PWGGA/PHOSTasks/CaloCellQA/macros
       )

install(FILES
                PHOS_TriggerQA/macros/AddTaskPHOSTriggerQA.C
        DESTINATION
                PWGGA/PHOSTasks/PHOS_TriggerQA/macros/
       )

install(FILES
                PHOS_PbPbQA/macros/AddTaskPHOSPbPb.C
        DESTINATION PWGGA/PHOSTasks/PHOS_PbPbQA/macros/
       )