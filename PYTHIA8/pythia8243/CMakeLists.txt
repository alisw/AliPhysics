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

# Definition of pythia8210. Please mind the naming convention:
#  - Module name: pythia8 (does not contain the minor version)
#    Use this when doing, i.e., set(LIBDEPS pythia8) from other modules.
#  - Library name: (on Unix) libpythia8243.so.
#    A symlink called libpythia8.so (pointing to the default minor version) is
#    also created for the convenience of the analyzers: i.e. you can just load
#    the default version with gSystem->Load("libpythia8").
set(MODULE pythia8)
add_definitions(-D_MODULE_="${MODULE}")
set(TARGETNAME pythia8243)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PYTHIA8/${TARGETNAME})

# Additional includes - alphabetical order except ROOT
include_directories(
  ${ROOT_INCLUDE_DIRS}
  ${AliRoot_SOURCE_DIR}/PYTHIA8/${TARGETNAME}/include
)

# Sources - alphabetical order
set(SRCS
	src/Analysis.cc
	src/Basics.cc
	src/BeamParticle.cc
	src/BeamRemnants.cc
	src/BeamShape.cc
	src/BoseEinstein.cc
	src/ColourReconnection.cc
	src/ColourTracing.cc
	src/DeuteronProduction.cc
	src/Event.cc
	src/FJcore.cc
	src/FragmentationFlavZpT.cc
	src/FragmentationSystems.cc
	src/GammaKinematics.cc
	src/HIUserHooks.cc
	src/HadronLevel.cc
	src/HadronScatter.cc
	src/HardDiffraction.cc
	src/HeavyIons.cc
	src/HelicityBasics.cc
	src/HelicityMatrixElements.cc
	src/HiddenValleyFragmentation.cc
	src/History.cc
	src/Info.cc
	src/JunctionSplitting.cc
	src/LHEF3.cc
	src/LesHouches.cc
	src/Merging.cc
	src/MergingHooks.cc
	src/MiniStringFragmentation.cc
	src/MultipartonInteractions.cc
	src/ParticleData.cc
	src/ParticleDecays.cc
	src/PartonDistributions.cc
	src/PartonLevel.cc
	src/PartonSystems.cc
	src/PartonVertex.cc
	src/PhaseSpace.cc
	src/ProcessContainer.cc
	src/ProcessLevel.cc
	src/Pythia.cc
	src/PythiaStdlib.cc
	src/RHadrons.cc
	src/ResonanceDecays.cc
	src/ResonanceWidths.cc
	src/ResonanceWidthsDM.cc
	src/Ropewalk.cc
	src/SLHAinterface.cc
	src/Settings.cc
	src/SigmaCompositeness.cc
	src/SigmaDM.cc
	src/SigmaEW.cc
	src/SigmaExtraDim.cc
	src/SigmaGeneric.cc
	src/SigmaHiggs.cc
	src/SigmaLeftRightSym.cc
	src/SigmaLeptoquark.cc
	src/SigmaNewGaugeBosons.cc
	src/SigmaOnia.cc
	src/SigmaProcess.cc
	src/SigmaQCD.cc
	src/SigmaSUSY.cc
	src/SigmaTotal.cc
	src/SimpleSpaceShower.cc
	src/SimpleTimeShower.cc
	src/SimpleWeakShowerMEs.cc
	src/StandardModel.cc
	src/Streams.cc
	src/StringFragmentation.cc
	src/StringLength.cc
	src/SusyCouplings.cc
	src/SusyLesHouches.cc
	src/SusyResonanceWidths.cc
	src/SusyWidthFunctions.cc
	src/TauDecays.cc
	src/UserHooks.cc
  )

# Note that we do not generate dictionary and .rootmap for Pythia.

# Add a library to the project using the specified source files
add_library_tested(${TARGETNAME} SHARED ${SRCS})
add_library(${MODULE} ALIAS ${TARGETNAME})

# Dependencies
set(LIBDEPS lhapdf)

# Additional compilation flags
set_target_properties(${TARGETNAME} PROPERTIES COMPILE_FLAGS "")
target_link_libraries(${TARGETNAME})

# Linking the library
target_link_libraries(${TARGETNAME} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${TARGETNAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${TARGETNAME} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Link this minor version to the default "pythia8". Respect library prefix and extension that
# depend on the destination platform
install(CODE "
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_SHARED_MODULE_PREFIX}${TARGETNAME}${CMAKE_SHARED_MODULE_SUFFIX}
    ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_MODULE_PREFIX}${MODULE}${CMAKE_SHARED_MODULE_SUFFIX}
  )
")

# Data needed to run Pythia. Install it with the explicit version name, and link it to "pythia8"
install(DIRECTORY share/Pythia8/xmldoc DESTINATION PYTHIA8/${TARGETNAME})
install(CODE "
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${TARGETNAME}
    ${CMAKE_INSTALL_PREFIX}/PYTHIA8/pythia8
  )
")
install(DIRECTORY include/Pythia8 DESTINATION include)
install(DIRECTORY include/Pythia8Plugins DESTINATION include)
