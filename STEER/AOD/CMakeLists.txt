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
set (MODULE AOD)

# Module include folders
include_directories(${CMAKE_SOURCE_DIR}/STEER/${MODULE})

# Additional include directories - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/STEER/ESD
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources - alphabetical order
set(SRCS
    AliAODBranchReplicator.cxx
    AliAODCaloCells.cxx
    AliAODCaloCluster.cxx
    AliAODCaloTrigger.cxx
    AliAODcascade.cxx
    AliAODCluster.cxx
    AliAODDiJet.cxx
    AliAODDimuon.cxx
    AliAODEvent.cxx
    AliAODExtension.cxx
    AliAODFmdCluster.cxx
    AliAODHandler.cxx
    AliAODHeader.cxx
    AliAODHMPIDrings.cxx
    AliAODInputHandler.cxx
    AliAODJet.cxx
    AliAODJetEventBackground.cxx
    AliAODMCHeader.cxx
    AliAODMCParticle.cxx
    AliAODPhoton.cxx
    AliAODPid.cxx
    AliAODpidUtil.cxx
    AliAODPmdCluster.cxx
    AliAODPWG4ParticleCorrelation.cxx
    AliAODPWG4Particle.cxx
    AliAODRecoDecay.cxx
    AliAODRedCov.cxx
    AliAODTagCreator.cxx
    AliAODTrack.cxx
    AliAODTracklets.cxx
    AliAODTrdTrack.cxx
    AliAODTrdTracklet.cxx
    AliAODTZERO.cxx
    AliAODv0.cxx
    AliAODVertex.cxx
    AliAODVZERO.cxx
    AliAODZDC.cxx
   )

# Sources from headers
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
