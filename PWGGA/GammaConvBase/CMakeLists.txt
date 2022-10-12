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
set(MODULE PWGGAGammaConvBase)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGGA/GammaConvBase)

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALtasks
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALtrigger
                    ${AliPhysics_SOURCE_DIR}/PWG/CaloTrackCorrBase
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools/yaml-cpp/include
                    ${AliPhysics_SOURCE_DIR}/PWG/TRD
                    ${AliPhysics_SOURCE_DIR}/TENDER/Tender
                    ${AliPhysics_SOURCE_DIR}/TENDER/TenderSupplies
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/ClusterSelection
                    ${AliPhysics_SOURCE_DIR}/PWGGA/PHOSTasks/PHOS_PbPb
                   )

# the KF particle library support if required
set(USEKFLIBRARY FALSE)
if(USEKFLIBRARY) 
    get_target_property(KFPARTICLE_INCLUDE_DIR KFParticle::KFParticle INTERFACE_INCLUDE_DIRECTORIES)
    include_directories(${KFPARTICLE_INCLUDE_DIR})
    add_definitions("-DWITH_KFPARTICLE")
endif(USEKFLIBRARY)

# Sources in alphabetical order
set(SRCS
    AliAODConversionMother.cxx
    AliAODConversionParticle.cxx
    AliAODConversionPhoton.cxx
    AliCaloPhotonCuts.cxx
    AliCaloSigmaCuts.cxx	
    AliCaloTrackMatcher.cxx
    AliConversionAODBGHandlerRP.cxx
    AliConversionCuts.cxx
    AliConversionMesonCuts.cxx
    AliConversionPhotonBase.cxx
    AliConversionPhotonCuts.cxx
    AliConversionSelection.cxx
    AliConversionTrackCuts.cxx
    AliConvEventCuts.cxx
    AliDalitzElectronCuts.cxx
    AliDalitzElectronSelector.cxx
    AliKFConversionMother.cxx
    AliKFConversionPhoton.cxx
    AliPhotonIsolation.cxx
    AliV0ReaderV1.cxx
    AliDalitzAODESD.cxx
    AliDalitzAODESDMC.cxx
    AliDalitzEventMC.cxx
    AliDalitzData.cxx
    AliGAKFParticle.cxx
    AliGAKFVertex.cxx
    AliResponseMatrixHelper.cxx
    AliGammaConvEventMixing.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS})
# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
if(USEKFLIBRARY)
    generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "-DPWGGAUSEKFPARTICLE")
else(USEKFLIBRARY)
    generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")
endif(USEKFLIBRARY)

set(ROOT_DEPENDENCIES Core EG GenVector Geom Gpad Hist MathCore Matrix Net Physics RIO Tree)
set(ALIROOT_DEPENDENCIES ANALYSIS ANALYSISalice AOD PWGEMCALtasks PWGEMCALbase PWGEMCALtrigger)
set(ALIPHYSICS_DEPENDENCIES EMCALbase PWGCaloTrackCorrBase PWGCaloTrackCorrBase OADB)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ALIPHYSICS_DEPENDENCIES} ${ROOT_DEPENDENCIES})
if(USEKFLIBRARY)
    get_target_property(KFPARTICLE_LIBRARY KFParticle::KFParticle IMPORTED_LOCATION)
    set(LIBDEPS ${LIBDEPS} ${KFPARTICLE_LIBRARY})
endif(USEKFLIBRARY)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)

# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${LIBDEPS})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs} INTERFACE ${KFPARTICLE_INCLUDE_DIR})

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
if(USEKFLIBRARY)
    target_compile_definitions(${MODULE} PUBLIC "PWGGAUSEKFPARTICLE")
endif(USEKFLIBRARY)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
