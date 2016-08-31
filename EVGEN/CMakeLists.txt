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
set(MODULE EVGEN)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/FASTSIM
                    ${AliRoot_SOURCE_DIR}/LHAPDF/base
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TEvtGen
                    ${AliRoot_SOURCE_DIR}/TEvtGen/THepMCParser
                    ${AliRoot_SOURCE_DIR}/PYTHIA6/AliPythia6
                   )

# Sources
set(SRCS
    AliDecayerExodus.cxx
    AliDecayerPolarized.cxx
    AliDimuCombinator.cxx
    AliGenAfterBurnerFlow.cxx
    AliGenBeamGas.cxx
    AliGenBox.cxx
    AliGenCocktailAfterBurner.cxx
    AliGenCocktail.cxx
    AliGenCocktailEntry.cxx
    AliGenCorrHF.cxx
    AliGenCosmicsParam.cxx
    AliGenDeuteron.cxx
    AliGenDoubleScan.cxx
    AliGenEposReader.cxx
    AliGenExtFile.cxx
    AliGenExtExec.cxx
    AliGenFixed.cxx
    AliGenFunction.cxx
    AliGenGeVSim.cxx
    AliGenGSIlib.cxx
    AliGenHalo.cxx
    AliGenHaloProtvino.cxx
    AliGenHBTosl.cxx
    AliGenHIJINGparaBa.cxx
    AliGenHIJINGpara.cxx
    AliGenHMPIDlib.cxx
    AliGenITSULib.cxx
    AliGenKrypton.cxx
    AliGenLcLib.cxx
    AliGenLib.cxx
    AliGenLightNuclei.cxx
    AliGenMC.cxx
    AliGenMUONCocktail.cxx
    AliGenMUONCocktailpp.cxx
    AliGenMUONlib.cxx
    AliGenMUONLMR.cxx
    AliGenPairFlat.cxx
    AliGenParam.cxx
    AliGenPHOSlib.cxx
    AliGenPileup.cxx
    AliGenPMDlib.cxx
    AliGenPromptPhotons.cxx
    AliGenReaderCwn.cxx
    AliGenReader.cxx
    AliGenReaderEcalHijing.cxx
    AliGenReaderEcalJets.cxx
    AliGenReaderEMD.cxx
    AliGenReaderHepMC.cxx
    AliGenReaderSL.cxx
    AliGenReaderTreeK.cxx
    AliGenScan.cxx
    AliGenSlowNucleons.cxx
    AliGenSTRANGElib.cxx
    AliGenThermalPhotons.cxx
    AliGenThetaSlice.cxx
    AliGenTHnSparse.cxx
    AliGenTunedOnPbPb.cxx
    AliGeVSimParticle.cxx
    AliOmegaDalitz.cxx
    AliSlowNucleonModel.cxx
    AliSlowNucleonModelExp.cxx
    AliGenPerformance.cxx
    AliGenReadersEMD.cxx
    AliGenToy.cxx
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Additional header files
set(HDRS ${HDRS} AliDecayer.h)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD TEvtGen FASTSIM THepMCParser)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

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
install(FILES ${HDRS} DESTINATION include)
install(DIRECTORY dataCorrHF DESTINATION EVGEN)

message(STATUS "${MODULE} enabled")
