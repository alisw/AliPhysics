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
set(MODULE PWGHFtreeHF)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGHF/treeHF)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
		    ${AliPhysics_SOURCE_DIR}/OADB
		    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
		    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Base
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Tasks
                    ${AliPhysics_SOURCE_DIR}/PWG/JETFW
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/TRD
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
		    ${AliPhysics_SOURCE_DIR}/PWGJE/EMCALJetTasks
  )

# Sources - alphabetical order
set(SRCS

  AliAnalysisTaskSEHFTreeCreator.cxx
  AliHFJet.cxx
  AliHFTreeHandler.cxx
  AliHFTreeHandlerD0toKpi.cxx
  AliHFTreeHandlerDplustoKpipi.cxx
  AliHFTreeHandlerDstoKKpi.cxx
  AliHFTreeHandlerLctopKpi.cxx
  AliHFTreeHandlerBplustoD0pi.cxx
  AliHFTreeHandlerBstoDspi.cxx
  AliHFTreeHandlerDstartoKpipi.cxx
  AliHFTreeHandlerLc2V0bachelor.cxx
  AliHFTreeHandlerLbtoLcpi.cxx
  AliHFTreeHandlerInclusiveJet.cxx
  AliJetTreeHandler.cxx
  AliParticleTreeHandler.cxx
  AliTrackletTreeHandler.cxx
  )

if(FASTJET_FOUND)
    set(SRCS ${SRCS}
        AliHFJetFinder.cxx
	)
endif(FASTJET_FOUND)

if(FASTJET_FOUND)
    include_directories(${AliPhysics_SOURCE_DIR}/JETAN/JETAN)
    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
    link_directories(${FASTJET_LIBS_DIR})
    add_definitions(${FASTJET_DEFINITIONS})
    list(APPEND HDRS FJ_includes.h)
endif(FASTJET_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${FASTJET_ROOTDICT_OPTS}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice PWGflowTasks PWGTRD PWGPPevcharQn PWGPPevcharQnInterface PWGHFvertexingHF PWGEMCALbase PWGJETFW PWGJEEMCALJetTasks PWGTools)

if(FASTJET_FOUND)
  set(LIBDEPS
    ${LIBDEPS}
    fastjetcontribfragile
    fastjetplugins
    siscone_spherical
    siscone
    fastjettools
    fastjet
    )
endif(FASTJET_FOUND)

generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependencies
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

install(DIRECTORY
            macros
        DESTINATION PWGHF/treeHF)

