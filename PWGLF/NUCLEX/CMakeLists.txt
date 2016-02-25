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
set(MODULE PWGLFnuclex)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/NUCLEX
                    ${AliPhysics_SOURCE_DIR}/PWGLF/NUCLEX/Utils/RecoDecay
                    ${AliPhysics_SOURCE_DIR}/PWGLF/NUCLEX/Utils/NuclexFilter
)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${ALIROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/Base
  )
set(SRCS
  Exotica/Hdibaryon/AliAnalysisTaskHdibaryonLPpi.cxx
  Exotica/LambdaN/AliAnalysisTaskLambdaNAOD.cxx
  Exotica/NOmega/AliAnalysisTaskNOmegaLPK.cxx
  Exotica/NOmega/AliAnalysisTaskOmegaOmegaOX.cxx
  Hypernuclei/Hyp2Body/AliAnalysisTaskHelium3Pi.cxx
  Hypernuclei/Hyp2Body/AliAnalysisTaskHelium3PiAOD.cxx
  Hypernuclei/Hyp2Body/AliAnalysisTaskHelium3PiMC.cxx
  Hypernuclei/Hyp3Body/AliAnalysisTaskHypertriton3.cxx
  Hypernuclei/Hyp3Body/AliAnalysisTaskHypertriton3AOD.cxx
  Nuclei/DeltaMasses/AliAnalysisMCNuclMult.cxx
  Nuclei/DeltaMasses/AliAnalysisNucleiMass.cxx
  Nuclei/DeltaMasses/AliAnalysisNuclMult.cxx
  Nuclei/DeuteronpA/AliAnalysisDeuteronpA.cxx
  Nuclei/DeuteronpA/AliAnalysisDeuteronTree.cxx
  Nuclei/He4/AliAnalysisTaskAntiHe4.cxx
  Nuclei/NucleiPbPb/AliAnalysisTaskNucleiYield.cxx
  Nuclei/Nucleipp/AliAnalysisTaskB2.cxx
  Nuclei/Nucleipp/AliAnalysisTaskB2AOD.cxx
  Nuclei/Nucleipp/AliLnAODtrackCuts.cxx
  Nuclei/Nucleipp/AliLnHistoMap.cxx
  Nuclei/Nucleipp/AliLnID.cxx
  Nuclei/v2/AliAnalysisTaskNucleiv2.cxx
  Nuclei/v2/AliAnalysisTaskNucleiv2SP.cxx
  Utils/NuclexFilter/AliAODMCNuclExReplicator.cxx
  Utils/NuclexFilter/AliAODNuclExReplicator.cxx
  Utils/NuclexFilter/AliAnalysisTaskESDNuclExFilter.cxx
  Utils/NuclexFilter/AliAnalysisTaskESDNuclExFilterMC.cxx
  Utils/NuclexFilter/AliAnalysisTaskReadNuclexAOD.cxx
  Utils/RecoDecay/AliAODRecoDecayLF.cxx
  Utils/RecoDecay/AliAODRecoDecayLF2Prong.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice CORRFW OADB PWGUDbase PWGHFvertexingHF)

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

# install macros
install(DIRECTORY . DESTINATION PWGLF/NUCLEX FILES_MATCHING PATTERN "*.C")

