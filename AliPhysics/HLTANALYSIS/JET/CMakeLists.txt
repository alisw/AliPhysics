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
set(MODULE AliHLTJET)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/HLTANALYSIS/JET
                    ${AliPhysics_SOURCE_DIR}/HLTANALYSIS/JET/analysis
                    ${AliPhysics_SOURCE_DIR}/HLTANALYSIS/JET/cone
                    ${AliPhysics_SOURCE_DIR}/HLTANALYSIS/JET/fastjet
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/JETAN/JETAN
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTJETAgent.cxx
    AliHLTJETBase.cxx
    AliHLTJETJetCuts.cxx
    AliHLTJETReader.cxx
    AliHLTJETReaderHeader.cxx
    AliHLTJETTrackCuts.cxx
    analysis/AliHLTJETAnalysisBase.cxx
    analysis/AliHLTJETAnalysisComponent.cxx
    analysis/AliHLTJETAnalysisJets.cxx
    analysis/AliHLTJETAnalysisMerge.cxx
    cone/AliHLTJETConeEtaPhiCell.cxx
    cone/AliHLTJETConeFinder.cxx
    cone/AliHLTJETConeGrid.cxx
    cone/AliHLTJETConeHeader.cxx
    cone/AliHLTJETConeJetCandidate.cxx
    cone/AliHLTJETConeJetComponent.cxx
    cone/AliHLTJETConeSeedCuts.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Check for FastJet installation to include proper headers
if(FASTJET_FOUND)
    set(HDRS ${HDRS}
        fastjet/AliHLTJETFastJetHeader.h
        fastjet/AliHLTJETFastJetFinder.h
        fastjet/AliHLTJETFastJetComponent.h
    )

    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
    link_directories(${FASTJET_LIBS_DIR})
    add_definitions(${FASTJET_DEFINITIONS})
endif(FASTJET_FOUND)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${FASTJET_ROOTDICT_OPTS}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase AOD ESD ANALYSISalice HLTbase AliHLTUtil JETAN)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")
target_link_libraries(${MODULE} STEERBase AOD ESD ANALYSISalice HLTbase AliHLTUtil JETAN)

if(FASTJET_FOUND)
    target_link_libraries(${MODULE} ${FASTJET_LIBS})
endif(FASTJET_FOUND)

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# Status message
message(STATUS "HLTJET enabled")
