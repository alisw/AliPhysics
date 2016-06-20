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

set(MODULE ConverterCollada)
# Module include folder
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/EVE/converters/${MODULE}
)

set(SRCS
    AliColladaBuffer.cxx
)

string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${MODULE_COMPILE_FLAGS}")

set(LIBDEPS EveBase Ged Gpad Graf Graf3d Gui Net RGL Rint RIO XMLIO)

generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")
# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${LIBDEPS})

set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")
# Installation
install(TARGETS ${MODULE}
ARCHIVE DESTINATION lib
LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
