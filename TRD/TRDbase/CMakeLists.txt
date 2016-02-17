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
set(MODULE TRDbase)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TRD/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliTRDalignment.cxx
    AliTRDarrayADC.cxx
    AliTRDarrayDictionary.cxx
    AliTRDarraySignal.cxx
    AliTRDCalChamberStatus.cxx
    AliTRDCalDCS.cxx
    AliTRDCalDCSFEE.cxx
    AliTRDCalDCSFEEv2.cxx
    AliTRDCalDCSGTUBoardInfo.cxx
    AliTRDCalDCSGTUCtpOpc.cxx
    AliTRDCalDCSGTU.cxx
    AliTRDCalDCSGTUSegment.cxx
    AliTRDCalDCSGTUTgu.cxx
    AliTRDCalDCSGTUTmu.cxx
    AliTRDCalDCSPTR.cxx
    AliTRDCalDCSv2.cxx
    AliTRDCalDet.cxx
    AliTRDCalibChamberStatus.cxx
    AliTRDcalibDB.cxx
    AliTRDCalibPadStatus.cxx
    AliTRDCalibraExbAltFit.cxx
    AliTRDCalibraFillHisto.cxx
    AliTRDCalibraFit.cxx
    AliTRDCalibraMode.cxx
    AliTRDCalibraVdriftLinearFit.cxx
    AliTRDCalibraVector.cxx
    AliTRDCalibViewer.cxx
    AliTRDCalibViewerGUI.cxx
    AliTRDCalMonitoring.cxx
    AliTRDCalOnlineGainTable.cxx
    AliTRDCalOnlineGainTableMCM.cxx
    AliTRDCalOnlineGainTableROC.cxx
    AliTRDCalPad.cxx
    AliTRDCalPadStatus.cxx
    AliTRDCalPID.cxx
    AliTRDCalROC.cxx
    AliTRDCalSingleChamberStatus.cxx
    AliTRDCalTrapConfig.cxx
    AliTRDCalTrkAttach.cxx
    AliTRDchamberTimeBin.cxx
    AliTRDCommonParam.cxx
    AliTRDcluster.cxx
    AliTRDclusterizer.cxx
    AliTRDDataDCS.cxx
    AliTRDdEdxBaseUtils.cxx
    AliTRDdEdxCalibHistArray.cxx
    AliTRDdEdxCalibUtils.cxx
    AliTRDdEdxReconUtils.cxx
    AliTRDdigit.cxx
    AliTRDdigitsManager.cxx
    AliTRDdigitsParam.cxx
    AliTRDEntriesInfo.cxx
    AliTRDfeeParam.cxx
    AliTRDgeometry.cxx
    AliTRDgtuParam.cxx
    AliTRDgtuSim.cxx
    AliTRDgtuTMU.cxx
    AliTRDhit.cxx
    AliTRDltuParam.cxx
    AliTRDmcmSim.cxx
    AliTRDonlineTrackMatching.cxx
    AliTRDpadPlane.cxx
    AliTRDPhInfo.cxx
    AliTRDpidUtil.cxx
    AliTRDPreprocessor.cxx
    AliTRDPrfInfo.cxx
    AliTRDQAChecker.cxx
    AliTRDrawData.cxx
    AliTRDrawStream.cxx
    AliTRDrecoParam.cxx
    AliTRDReconstructor.cxx
    AliTRDSaxHandler.cxx
    AliTRDseedV1.cxx
    AliTRDSensorArray.cxx
    AliTRDSensor.cxx
    AliTRDSignalIndex.cxx
    AliTRDSimParam.cxx
    AliTRDtrackGTU.cxx
    AliTRDtrackerDebug.cxx
    AliTRDtrackerV1.cxx
    AliTRDtrackingChamber.cxx
    AliTRDtrackingSector.cxx
    AliTRDtrackletBase.cxx
    AliTRDtrackletGTU.cxx
    AliTRDtrackletOflHelper.cxx
    AliTRDtrackletMCM.cxx
    AliTRDtrackletWord.cxx
    AliTRDtrackOnline.cxx
    AliTRDtrackV1.cxx
    AliTRDtransform.cxx
    AliTRDtrapConfig.cxx
    AliTRDtrapConfigHandler.cxx
    AliTRDUshortInfo.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Setting compile and link flags in case of date
if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS} ${DATE_DYNMONLIBS}")
endif(DATE_FOUND)

# Check if CRS support was enabled
if(TRD_RAW_CRC)
    message(STATUS "CRC support for AliTRDrawStream enabled")
    set(MODULE_COMPILE_FLAGS "-DTRD_RAW_CRC ${MODULE_COMPILE_FLAGS}")
endif(TRD_RAW_CRC)

# Check if RAW debug enabled
if(TRD_RAW_DEBUG)
    message(STATUS "Debug support for AliTRDrawStream enabled")
    set(MODULE_COMPILE_FLAGS "-DTRD_RAW_DEBUG ${MODULE_COMPILE_FLAGS}")
endif(TRD_RAW_DEBUG)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${MODULE_COMPILE_FLAGS}")

set(ROOT_DEPENDENCIES Core GenVector Geom Gpad Graf Gui Hist MathCore Matrix Minuit Net RIO Tree VMC XMLParser)
set(ALIROOT_DEPENDENCIES STEERBase STEER ESD CDB RAWDatabase RAWDatarec)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Follow headers dependencies
add_dependencies(${MODULE}-object ${ALIROOT_DEPENDENCIES})
# Add a library to the project using the object
add_library(${MODULE} SHARED $<TARGET_OBJECTS:TRDbase-object>)

# Linking
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)


# Setting the correct headers for the object as gathered from the dependencies
# Additional compilation flags
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Setting the correct headers for the object as gathered from the dependencies
# This has to be added after all the compile definitions
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# Static version if DA enabled
if(ALIROOT_STATIC)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:TRDbase-object>)
    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})

    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")
    target_link_libraries(${MODULE}-static ${STATIC_ALIROOT_DEPENDENCIES} Root RootExtra)
    
    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})
    
    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endif(ALIROOT_STATIC)
