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
set(MODULE STEERBase)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/STEER/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_BINARY_DIR}/version
                   )

# Sources - alphabetical order
set(SRCS
    AliCentrality.cxx
    AliCheb2DStack.cxx
    AliCheb2DStackF.cxx
    AliCheb2DStackS.cxx
    AliCheb3DCalc.cxx
    AliCheb3D.cxx
    AliCodeTimer.cxx
    AliCollisionGeometry.cxx
    AliCounterCollection.cxx
    AliDAQ.cxx
    AliDataFile.cxx
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
    AliFileUtilities.cxx
    AliGenCocktailEventHeader.cxx
    AliGenDPMjetEventHeader.cxx
    AliGenEpos3EventHeader.cxx
    AliGenEposEventHeader.cxx
    AliGenEventHeader.cxx
    AliGenEventHeaderTunedPbPb.cxx
    AliGenGeVSimEventHeader.cxx
    AliGenHepMCEventHeader.cxx
    AliGenHerwigEventHeader.cxx
    AliGenHijingEventHeader.cxx
    AliGenPythiaEventHeader.cxx
    AliGenToyEventHeader.cxx
    AliHeader.cxx
    AliHLTBaseVGRPAccess.cxx
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
    AliVAD.cxx
    AliVAODHeader.cxx
    AliVCaloCells.cxx
    AliVCaloTrigger.cxx
    AliVCluster.cxx
    AliVCuts.cxx
    AliVersion.cxx
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
    AliVTPCseed.cxx
    AliVTrack.cxx
    AliVTrdTrack.cxx
    AliVTrdTracklet.cxx
    AliVVertex.cxx
    AliVVZERO.cxx
    AliVVZEROfriend.cxx
    AliVZDC.cxx
    TTreeStream.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

set(HDRS ${HDRS} 
    AliSharedPtr.h
    )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(STEERBase-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library(STEERBase SHARED $<TARGET_OBJECTS:STEERBase-object>)

# Generate the ROOT map
# Dependecies
set(LIBDEPS Core EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Net Physics RIO Tree VMC)
generate_rootmap("STEERBase" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Additional compilation flags for the object
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "")
# Linking the library, not the object
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(STEERBase PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Static version if DA enabled
if(ALIROOT_STATIC)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:STEERBase-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})
    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")
    target_link_libraries(${MODULE}-static Root RootExtra)

    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)



#TEST programs
add_executable(TestAliFileUtilities TestAliFileUtilities.cxx)
target_link_libraries(TestAliFileUtilities STEERBase ${ROOT_LIBRARIES})
