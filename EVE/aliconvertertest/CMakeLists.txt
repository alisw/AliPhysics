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

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/EVE/aliconvertertest
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/EVE/Converter
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )
# alistorage executable
add_executable(aliconvertertest aliconvertertest.cxx)
target_link_libraries(aliconvertertest Converter STEERBase STEER ESD)

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

# Setting target compile and linking flags
set_target_properties(aliconvertertest PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(aliconvertertest PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

if(PYTHIA6)
    target_link_libraries(aliconvertertest EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS aliconvertertest RUNTIME DESTINATION bin)
