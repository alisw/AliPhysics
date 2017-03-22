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
set(MODULE StarLight)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/STARLIGHT/starlight/include
                    ${CMAKE_CURRENT_BINARY_DIR}
)

# Additional include folders in alphabetical order
include_directories(
                   )

# Sources in alphabetical order
set(SRCS
    src/beambeamsystem.cpp
    src/beam.cpp
    src/bessel.cpp
    src/eventchannel.cpp
    src/eventfilewriter.cpp
    src/filewriter.cpp
    src/gammaaluminosity.cpp
    src/gammaavm.cpp
    src/gammagammaleptonpair.cpp
    src/gammagammasingle.cpp
    src/incoherentPhotonNucleusLuminosity.cpp
    src/incoherentVMCrossSection.cpp
    src/inputParameters.cpp
    src/inputParser.cpp
    src/lorentzvector.cpp
    src/narrowResonanceCrossSection.cpp
    src/nBodyPhaseSpaceGen.cpp
    src/nucleus.cpp
    src/photonNucleusCrossSection.cpp
    src/randomgenerator.cpp
    src/readinluminosity.cpp
    src/starlight.cpp
    src/starlightparticlecodes.cpp
    src/starlightparticle.cpp
    src/starlightStandalone.cpp
    src/twophotonluminosity.cpp
    src/upcevent.cpp
    src/vector3.cpp
    src/wideResonanceCrossSection.cpp
   )

set(Starlight_VERSION_MAJOR 9999)
set(Starlight_VERSION_MINOR 1)
set(Starlight_VERSION_MINOR_MINOR 0)

configure_file(${AliRoot_SOURCE_DIR}/STARLIGHT/starlight/starlightconfig.h.in
               starlightconfig.h
              )
# Headers in alphabetical order
set(HDRS
    include/beambeamsystem.h
    include/beam.h
    include/bessel.h
    include/eventchannel.h
    include/eventfilewriter.h
    include/filewriter.h
    include/gammaaluminosity.h
    include/gammaavm.h
    include/gammagammaleptonpair.h
    include/gammagammasingle.h
    include/incoherentPhotonNucleusLuminosity.h
    include/incoherentVMCrossSection.h
    include/inputParameters.h
    include/inputParser.h
    include/lorentzvector.h
    include/narrowResonanceCrossSection.h
    include/nBodyPhaseSpaceGen.h
    include/nucleus.h
    include/photonNucleusCrossSection.h
    include/pythiadecayer.h
    include/PythiaStarlight.h
    include/randomgenerator.h
    include/readinluminosity.h
    include/reportingUtils.h
    include/spectrum.h
    include/spectrumprotonnucleus.h
    include/starlightconstants.h
    include/starlightdpmjet.h
    include/starlight.h
    include/starlightlimits.h
    include/starlightparticlecodes.h
    include/starlightparticle.h
    include/starlightpythia.h
    include/starlightStandalone.h
    include/twophotonluminosity.h
    include/upcevent.h
    include/vector3.h
    include/wideResonanceCrossSection.h

   )

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS})

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

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/starlightconfig.h DESTINATION include)
install(FILES ${HDRS} DESTINATION include)
