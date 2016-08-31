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
set(MODULE MONITORzmq)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/${MODULE}
                    ${AliRoot_SOURCE_DIR}/MONITOR/alistoragemanager
                    ${AliRoot_SOURCE_DIR}/MONITOR/alionlinereco
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/MONITOR
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCrec
                    ${ZEROMQ_INCLUDE_DIR}
                   )

# Sources in alphabetical order
set(SRCS
    AliCommunicationThread.cxx
    AliDIMListenerThread.cxx
    AliEventsCollectorThread.cxx
    AliOnlineReconstruction.cxx
    AliStorageAdministratorPanel.cxx
    AliStorageAdministratorPanelListEvents.cxx
    AliStorageAdministratorPanelMarkEvent.cxx
    AliStorageAdministratorPanelSetStorageParams.cxx
    AliStorageClientThread.cxx
    AliStorageDatabase.cxx
    AliZMQManager.cxx
    AliStorageServerThread.cxx
   )

if(DIM_FOUND)
    set(SRCS ${SRCS}
        ${AliRoot_SOURCE_DIR}/MONITOR/MONITOR/SORNotifier.cxx
       )
    include_directories(${DIM_INCLUDE_DIR})
    add_definitions(-DALI_DIM)
endif()

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS} ${DATE_DYNMONLIBS}")
endif(DATE_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${MODULE_COMPILE_FLAGS}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS TPCrec STEERBase STEER ESD CDB MONITOR ${ZEROMQ_LIBRARIES})
if(DIM_FOUND)
    set(LIBDEPS ${LIBDEPS} ${DIM_LIBRARIES})
endif()
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Additional compilation and linking flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
