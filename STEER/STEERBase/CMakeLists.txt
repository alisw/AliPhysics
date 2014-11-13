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
set (MODULE STEERBase)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/STEER/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS})

# Sources - alphabetical order
set(SRCS
    AliCentrality.cxx
    AliCheb3DCalc.cxx
    AliCheb3D.cxx
    AliCodeTimer.cxx
    AliCollisionGeometry.cxx
    AliCounterCollection.cxx
    AliDAQ.cxx
    AliDetectorEventHeader.cxx
    AliDetectorPID.cxx
    AliDetectorTagCuts.cxx
    AliDetectorTag.cxx
    AliDummyHandler.cxx
    AliEMCALPIDResponse.cxx
    AliEventplane.cxx
    AliEventTagCuts.cxx
    AliEventTag.cxx
    AliExternalTrackParam.cxx
    AliFileTag.cxx
    AliGenCocktailEventHeader.cxx
    AliGenDPMjetEventHeader.cxx
    AliGenEposEventHeader.cxx
    AliGenEventHeader.cxx
    AliGenEventHeaderTunedPbPb.cxx
    AliGenGeVSimEventHeader.cxx
    AliGenHepMCEventHeader.cxx
    AliGenHerwigEventHeader.cxx
    AliGenHijingEventHeader.cxx
    AliGenPythiaEventHeader.cxx
    AliHeader.cxx
    AliHMPIDPIDParams.cxx
    AliHMPIDPIDResponse.cxx
    AliInputEventHandler.cxx
    AliITSPidParams.cxx
    AliITSPIDResponse.cxx
    AliLHCTagCuts.cxx
    AliLHCTag.cxx
    AliLog.cxx
    AliMagF.cxx
    AliMagWrapCheb.cxx
    AliMathBase.cxx
    AliMCEvent.cxx
    AliMCEventHandler.cxx
    AliMCParticle.cxx
    AliMCVertex.cxx
    AliMergeableCollection.cxx
    AliMixedEvent.cxx
    AliNeutralTrackParam.cxx
    AliOADBContainer.cxx
    AliPDG.cxx
    AliPIDCombined.cxx
    AliPID.cxx
    AliPIDResponse.cxx
    AliPIDValues.cxx
    AliProdInfo.cxx
    AliQA.cxx
    AliRefArray.cxx
    AliRunTagCuts.cxx
    AliRunTag.cxx
    AliStack.cxx
    AliSysInfo.cxx
    AliTagCreator.cxx
    AliTimeStamp.cxx
    AliTOFHeader.cxx
    AliTOFPIDParams.cxx
    AliTOFPIDResponse.cxx
    AliTOFTriggerMask.cxx
    AliTPCdEdxInfo.cxx
    AliTPCPIDResponse.cxx
    AliTrackReference.cxx
    AliTRDdEdxParams.cxx
    AliTRDNDFast.cxx
    AliTRDPIDParams.cxx
    AliTRDPIDReference.cxx
    AliTRDPIDResponse.cxx
    AliTRDPIDResponseObject.cxx
    AliTRDTKDInterpolator.cxx
    AliTriggerScalers.cxx
    AliTriggerScalersRecord.cxx
    AliVAODHeader.cxx
    AliVCaloCells.cxx
    AliVCaloTrigger.cxx
    AliVCluster.cxx
    AliVCuts.cxx
    AliVEvent.cxx
    AliVEventHandler.cxx
    AliVEventPool.cxx
    AliVfriendEvent.cxx
    AliVfriendTrack.cxx
    AliVHeader.cxx
    AliVMFT.cxx
    AliVMisc.cxx
    AliVMultiplicity.cxx
    AliVParticle.cxx
    AliVTOFcluster.cxx
    AliVTOFHit.cxx
    AliVTOFMatch.cxx
    AliVTrack.cxx
    AliVTrdTrack.cxx
    AliVTrdTracklet.cxx
    AliVVertex.cxx
    AliVVZERO.cxx
    AliVZDC.cxx
    TTreeStream.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS Core EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Net Physics RIO Tree VMC)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "-O -g")
target_link_libraries(${MODULE} Core EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Net Physics RIO Tree VMC)

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
