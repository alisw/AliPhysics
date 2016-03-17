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
set(MODULE PWGLFppvsmult)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/ppVsMult)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
  )

# Sources - alphabetical order
set(SRCS
  dNdEtaPP13/task/AliITSMultRecBg.cxx
  dNdEtaPP13/task/AliAnalysisTaskdNdEtapp13.cxx
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Extra files to include in the PARfile
#set(PAR_EXTRAFILES
#
#)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependencies
set(LIBDEPS ANALYSISalice PWGTools ITSbase)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
#add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "PiKaPr/TestAOD PiKaPr/TPCTOFfits XtAnalysis" "${PAR_EXTRAFILES}")

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
install(DIRECTORY dNdEtaPP13/macros DESTINATION  PWGLF/ppVsMult/dNdEtaPP13 FILES_MATCHING PATTERN "*.C")
install(DIRECTORY dNdEtaPP13/task DESTINATION  PWGLF/ppVsMult/dNdEtaPP13 FILES_MATCHING PATTERN "*.C")
#install(DIRECTORY ChargedHadrons/dNdPt/MeanPtRaw/macros DESTINATION
#  PWGLF/SPECTRA/ChargedHadrons/dNdPt/MeanPtRaw)
#install(DIRECTORY IdentifiedHighPt/train/macros DESTINATION  PWGLF/SPECTRA/IdentifiedHighPt/train)
#install(DIRECTORY Kinks/macros DESTINATION  PWGLF/SPECTRA/Kinks)
#install(DIRECTORY Nuclei/B2/macros DESTINATION PWGLF/SPECTRA/Nuclei/B2)
#install(DIRECTORY Nuclei/masses  DESTINATION PWGLF/SPECTRA/Nuclei)
#install(DIRECTORY PiKaPr/TestAOD DESTINATION PWGLF/SPECTRA/PiKaPr)
#install(DIRECTORY PiKaPr/ITSsa DESTINATION PWGLF/SPECTRA/PiKaPr)
#install(DIRECTORY PiKaPr/TPCTOFfits DESTINATION PWGLF/SPECTRA/PiKaPr)
#install(FILES PiKaPr/TPCTOF/AddTaskCombinedHadronSpectra.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOF)
#install(FILES PiKaPr/TPCTOFpA/AddTaskTPCTOFpA.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOFpA)
#install(DIRECTORY XtAnalysis/macros DESTINATION PWGLF/SPECTRA/XtAnalysis)
#install(DIRECTORY MultEvShape DESTINATION PWGLF/SPECTRA FILES_MATCHING PATTERN "*.C")
#install(FILES MultEvShape/priorsDist_data_LHC10d_newAliroot.root DESTINATION  PWGLF/SPECTRA/MultEvShape)
#install(FILES Spherocity/macros/AddTaskEventShape.C DESTINATION  PWGLF/SPECTRA/Spherocity/macros)
#install(DIRECTORY AntiprotonToProton/train DESTINATION PWGLF/SPECTRA/AntiprotonToProton)
