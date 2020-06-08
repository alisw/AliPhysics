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
set(MODULE PWGMMue)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGMM/UE)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGLF/SPECTRA/Spherocity
  )

# Sources - alphabetical order
set(SRCS
  UeSpectra/AliAnalysisTaskPPvsUE.cxx
  UeSpherocity/AliAnalysisTaskUeSpherocity.cxx
  UeSpherocity/AliAnalysisTaskGenUeSpherocity.cxx
  UeSpherocity/AliAnalysisTaskUeRSpherocity.cxx
  UePid/AliAnalysisTaskUePid.cxx
  UeMultDep/AliAnalysisTaskUeMultDep.cxx
  UeDphi/AliAnalysisTaskUeSpectraDphi.cxx
  UeSpectraRT/AliAnalysisTaskUeSpectraRT.cxx
  UeKNO/AliAnalysisTaskMcKnoUe.cxx
  UeKNO/AliAnalysisTaskMcKno.cxx
  UeKNO/AliAnalysisTaskDataKno.cxx
  UeKNO/AliAnalysisTaskGenMcKno.cxx
  UeKNO/AliAnalysisTaskGenUeNchTS.cxx
  UeMpi/AliAnalysisTaskMpiUE.cxx		
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
	generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

	set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES ANALYSISalice PWGTools)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
	generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library. Note the extra includes ("Tracks UserTasks")
	add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "Tracks UserTasks")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)

# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

	set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# Additional compilation and linking flags
	set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
	set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting compilation flags for the object
	set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
	set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
	install(TARGETS ${MODULE}
			ARCHIVE DESTINATION lib
			LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# install macros
	install(DIRECTORY UeSpectra DESTINATION  PWGMM/UE)
	install(DIRECTORY UeSpherocity DESTINATION  PWGMM/UE)
	install(DIRECTORY UePid DESTINATION  PWGMM/UE)
        install(DIRECTORY UeKNO DESTINATION  PWGMM/UE)
        install(DIRECTORY UeMpi DESTINATION  PWGMM/UE)
install(DIRECTORY UeMultDep DESTINATION  PWGMM/UE)
install(DIRECTORY UeSpectraRT DESTINATION  PWGMM/UE)
install(DIRECTORY UeDphi DESTINATION  PWGMM/UE)
install(DIRECTORY macros DESTINATION  PWGMM/UE)		
