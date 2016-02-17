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
set(MODULE EveDet)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/EVE/${MODULE})

# Additional include folders in alphabetical order
include_directories(
                    ${ROOT_INCLUDE_DIR}
		    ${AliRoot_SOURCE_DIR}/AD/ADbase
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALbase
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${AliRoot_SOURCE_DIR}/EVE/EveBase
                    ${AliRoot_SOURCE_DIR}/FMD/FMDbase
                    ${AliRoot_SOURCE_DIR}/FMD/FMDrec
                    ${AliRoot_SOURCE_DIR}/FMD/FMDsim
                    ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                    ${AliRoot_SOURCE_DIR}/ITS/ITSsim
  	            ${AliRoot_SOURCE_DIR}/ITSMFT/COMMONITSMFT/COMMONITSMFTBase
                    ${AliRoot_SOURCE_DIR}/ITSMFT/ITS/ITSUpgradeBase
                    ${AliRoot_SOURCE_DIR}/MONITOR/alistoragemanager
                    ${AliRoot_SOURCE_DIR}/MUON/MUONbase
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON/MUONgeometry
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/MUON/MUONrec
                    ${AliRoot_SOURCE_DIR}/MUON/MUONsim
                    ${AliRoot_SOURCE_DIR}/PMD/PMDbase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/T0/T0base
                    ${AliRoot_SOURCE_DIR}/TOF/TOFbase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                    ${AliRoot_SOURCE_DIR}/TRD/TRDbase
                    ${AliRoot_SOURCE_DIR}/TRD/TRDrec
                    ${AliRoot_SOURCE_DIR}/VZERO/VZEROrec
                   )

# Sources in alphabetical order
set(SRCS
    AliEveADModule.cxx
    AliEveEMCALData.cxx
    AliEveEMCALSModule.cxx
    AliEveEMCALSModuleData.cxx
    AliEveFMDLoader.cxx
    AliEveITSDigitsInfo.cxx
    AliEveITSModule.cxx
    AliEveITSModuleStepper.cxx
    AliEveITSModuleStepperEditor.cxx
    AliEveITSScaledModule.cxx
    AliEveITSScaledModuleEditor.cxx
    AliEveITSUModule.cxx
    AliEveListAnalyser.cxx
    AliEveListAnalyserEditor.cxx
    AliEveMUONChamber.cxx
    AliEveMUONChamberData.cxx
    AliEveMUONChamberEditor.cxx
    AliEveMUONChamberGL.cxx
    AliEveMUONData.cxx
    AliEveMUONTrack.cxx
    AliEvePMDModule.cxx
    AliEvePMDModuleEditor.cxx
    AliEveT0Module.cxx
    AliEveTOFDigitsInfo.cxx
    AliEveTOFDigitsInfoEditor.cxx
    AliEveTOFSector.cxx
    AliEveTOFSectorEditor.cxx
    AliEveTOFStrip.cxx
    AliEveTOFStripEditor.cxx
    AliEveTPCData.cxx
    AliEveTPCLoader.cxx
    AliEveTPCLoaderEditor.cxx
    AliEveTPCSector2D.cxx
    AliEveTPCSector2DEditor.cxx
    AliEveTPCSector2DGL.cxx
    AliEveTPCSector3D.cxx
    AliEveTPCSector3DEditor.cxx
    AliEveTPCSector3DGL.cxx
    AliEveTPCSectorData.cxx
    AliEveTPCSectorViz.cxx
    AliEveTPCSectorVizEditor.cxx
    AliEveTRDData.cxx
    AliEveTRDDataEditor.cxx
    AliEveTRDLoader.cxx
    AliEveTRDLoaderImp.cxx
    AliEveTRDLoaderManager.cxx
    AliEveTRDModule.cxx
    AliEveTRDModuleImp.cxx
    AliEveTRDTrackList.cxx
    AliEveTRDTrackListEditor.cxx
    AliEveVZEROModule.cxx
    AliEveVZEROModuleEditor.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD CDB RAWDatabase RAWDatarec ADbase EMCALbase EMCALUtils FMDbase FMDrec FMDsim ITSbase ITSsim ITSUpgradeBase ITSUpgradeSim MUONbase MUONgeometry MUONmapping MUONrec MUONsim PMDbase T0base TOFbase TPCbase TRDbase TRDrec VZERObase VZEROrec EveBase EG Eve Ged GenVector Gpad Gui Minuit RGL STAT MLP)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} STEERBase STEER ESD CDB RAWDatabase RAWDatarec ADbase EMCALbase EMCALUtils FMDbase FMDrec FMDsim ITSbase ITSsim ITSUpgradeBase ITSUpgradeSim MUONbase MUONgeometry MUONmapping MUONrec MUONsim PMDbase T0base TOFbase TPCbase TRDbase TRDrec VZERObase VZEROrec EveBase EG Eve Ged GenVector Gpad Gui Minuit RGL STAT MLP)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
