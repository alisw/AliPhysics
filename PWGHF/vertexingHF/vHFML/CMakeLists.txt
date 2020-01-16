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

if(ROOT_VERSION_MAJOR EQUAL 6)
# Module
set(MODULEHFML PWGHFvertexingHFML)
add_definitions(-D_MODULEHFML_="${MODULEHFML}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF/vHFML)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/ML
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF/charmFlow
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Base
)

# Sources - alphabetical order
set(SRCS 
    AliHFMLResponse.cxx
    AliHFMLResponseDplustoKpipi.cxx
    AliHFMLResponseDstoKKpi.cxx
    AliHFMLResponseLctoV0bachelor.cxx
    AliHFMLVarHandler.cxx
    AliHFMLVarHandlerDplustoKpipi.cxx
    AliHFMLVarHandlerDstoKKpi.cxx
    AliAnalysisTaskSEDplus.cxx
    AliAnalysisTaskSEDs.cxx
    AliAnalysisTaskSECharmHadronvn.cxx
    AliAnalysisTaskSEHFTreeCreatorApply.cxx
    AliHFTreeHandlerApply.cxx
    AliHFTreeHandlerApplyDstoKKpi.cxx
    AliHFTreeHandlerApplyLc2V0bachelor.cxx
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULEHFML}" "${MODULEHFML}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULEHFML} SHARED  ${SRCS} G__${MODULEHFML}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ML PWGHFvertexingHF)
generate_rootmap("${MODULEHFML}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULEHFML}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULEHFML} "${SRCS}" "${HDRS}" "${MODULEHFML}LinkDef.h" "${LIBDEPS}")

# Linking the library
target_link_libraries(${MODULEHFML} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULEHFML} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULEHFML} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULEHFML}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

endif()
