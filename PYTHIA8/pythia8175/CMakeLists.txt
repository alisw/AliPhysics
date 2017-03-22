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
set(MODULE pythia8175)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PYTHIA8/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/PYTHIA8/${MODULE}/include
		   )

# Sources - alphabetical order
set(SRCS
  src/Analysis.cxx 
  src/Basics.cxx 
  src/BeamParticle.cxx 
  src/BeamRemnants.cxx 
  src/BeamShape.cxx 
  src/BoseEinstein.cxx 
  src/Event.cxx 
  src/FragmentationFlavZpT.cxx 
  src/FragmentationSystems.cxx 
  src/HadronLevel.cxx 
  src/HadronScatter.cxx
  src/HelicityBasics.cxx
  src/HelicityMatrixElements.cxx
  src/HiddenValleyFragmentation.cxx
  src/History.cxx
  src/Info.cxx 
  src/LesHouches.cxx 
  src/Merging.cxx 
  src/MergingHooks.cxx 
  src/MiniStringFragmentation.cxx 
  src/MultipartonInteractions.cxx 
  src/ParticleData.cxx 
  src/ParticleDecays.cxx 
  src/PartonDistributions.cxx 
  src/PartonLevel.cxx 
  src/PartonSystems.cxx 
  src/PhaseSpace.cxx 
  src/ProcessContainer.cxx 
  src/ProcessLevel.cxx 
  src/Pythia.cxx 
  src/PythiaStdlib.cxx
  src/RHadrons.cxx 
  src/ResonanceDecays.cxx 
  src/ResonanceWidths.cxx 
  src/Settings.cxx 
  src/SigmaCompositeness.cxx 
  src/SigmaEW.cxx 
  src/SigmaExtraDim.cxx 
  src/SigmaGeneric.cxx 
  src/SigmaHiggs.cxx 
  src/SigmaLeftRightSym.cxx 
  src/SigmaLeptoquark.cxx 
  src/SigmaNewGaugeBosons.cxx 
  src/SigmaOnia.cxx 
  src/SigmaProcess.cxx 
  src/SigmaQCD.cxx 
  src/SigmaSUSY.cxx 
  src/SigmaTotal.cxx 
  src/SpaceShower.cxx 
  src/StandardModel.cxx 
  src/StringFragmentation.cxx 
  src/SusyCouplings.cxx 
  src/SusyLesHouches.cxx 
  src/TimeShower.cxx 
  src/TauDecays.cxx 	  	  
  src/UserHooks.cxx 
  src/SusyResonanceWidths.cxx 
  )

# Headers generated from sources list
# NOTE: this is non-standard and requires two replaces
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
string(REPLACE "src/" "include/" HDRS "${HDRS}")

# Additional headers, with no corresponding sources
list(APPEND HDRS "include/FastJet3.h"
                 "include/HepMCInterface.h"
                 "include/LHAFortran.h"
                 "include/LHAPDFInterface.h"
                 "include/PythiaComplex.h"
    )

# No dictionary

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS})

# Do not generate the ROOT map
# Dependencies
set(LIBDEPS)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")
target_link_libraries(${MODULE} )

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

# Installation of headers
# NOTE: they will all go in the include/pythia subdirectory, and NOT in include
install(FILES ${HDRS} DESTINATION "include/pythia")
