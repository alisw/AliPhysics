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

# Checks for a ZeroMQ installation. Enables ZeroMQ by default if found on the
# system and with the right version.
#
# Point to a custom ZeroMQ installation with -DZEROMQ=<path>: in that case, if
# ZeroMQ is not found or has not the right version, a fatal error is raised.
#
# Variables set:
#
#  - ZEROMQ_FOUND          True if ZeroMQ is found
#  - ZEROMQ_INCLUDE_DIR    Where to find ZeroMQ include directory
#  - ZEROMQ_LIBRARIES      List of libraries when using ZeroMQ
#  - ZEROMQ_VERSION        ZeroMQ version, major.minor.patches
#  - ZEROMQ_VERSION_MAJOR  Major component of version
#  - ZEROMQ_VERSION_MINOR  Minor component of version
#  - ZEROMQ_VERSION_MAJOR  Patches component of version

set(ZEROMQ_VERSION_MIN "4.0.0")

message(STATUS "Checking for ZeroMQ ${ZEROMQ}")

set(ZEROMQ_FOUND FALSE)

if(ZEROMQ)
    # ZeroMQ is installed in a custom place
    find_library(ZEROMQ_LIBRARIES NAMES zmq
                PATHS ${ZEROMQ}/lib
                NO_DEFAULT_PATH
                DOC "Path to libzmq)"
            )
    find_path(ZEROMQ_INCLUDE_DIR NAMES zmq.h zmq_utils.h
                PATHS ${ZEROMQ}/include
                NO_DEFAULT_PATH
                DOC "Path to ZeroMQ include header files."
            )       
else(ZEROMQ)
        # Check is the library is installed on the system
    find_library(ZEROMQ_LIBRARIES NAMES zmq
                DOC "Path to libzmq)"
            )

    find_path(ZEROMQ_INCLUDE_DIR NAMES zmq_utils.h
                DOC "Path to ZeroMQ include header files."
            )
endif(ZEROMQ)

mark_as_advanced(ZEROMQ_LIBRARIES ZEROMQ_INCLUDE_DIR)

# Did we find ZeroMQ? We need to parse its version.
if(ZEROMQ_INCLUDE_DIR)
  message(STATUS "ZeroMQ include path: ${ZEROMQ_INCLUDE_DIR}")
endif()
if(ZEROMQ_LIBRARIES)
  message(STATUS "ZeroMQ libraries: ${ZEROMQ_LIBRARIES}")
endif()
if(ZEROMQ_INCLUDE_DIR AND ZEROMQ_LIBRARIES)
  file(READ "${ZEROMQ_INCLUDE_DIR}/zmq.h" zmqh)

  string(REGEX MATCH "#define +ZMQ_VERSION_MAJOR +([0-9]+)" zmqv "${zmqh}")
  set(ZEROMQ_VERSION_MAJOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#define +ZMQ_VERSION_MINOR +([0-9]+)" zmqv "${zmqh}")
  set(ZEROMQ_VERSION_MINOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#define +ZMQ_VERSION_PATCH +([0-9]+)" zmqv "${zmqh}")
  set(ZEROMQ_VERSION_PATCH "${CMAKE_MATCH_1}")

  unset(zmqh)
  unset(zmqv)

  set(ZEROMQ_VERSION "${ZEROMQ_VERSION_MAJOR}.${ZEROMQ_VERSION_MINOR}.${ZEROMQ_VERSION_PATCH}")
  message(STATUS "ZeroMQ version: ${ZEROMQ_VERSION}")

  if(ZEROMQ_VERSION VERSION_GREATER "${ZEROMQ_VERSION_MIN}")
    # Version OK.
    set(ZEROMQ_FOUND TRUE)
    add_definitions(-DZMQ)
    message(STATUS "ZeroMQ version ${ZEROMQ_VERSION} (> ${ZEROMQ_VERSION_MIN}) found")
  elseif(ZEROMQ)
    # Version not OK and explicitly requested: fatal.
    message(FATAL_ERROR "ZeroMQ in ${ZEROMQ} has version ${ZEROMQ_VERSION} <= ${ZEROMQ_VERSION_MIN}")
  else()
    message(STATUS "ZeroMQ found but version ${ZEROMQ_VERSION} <= ${ZEROMQ_VERSION_MIN}. Disabling ZeroMQ support.")
  endif()
elseif(ZEROMQ)
  # ZeroMQ not found and explicitly requested: fatal.
  message(FATAL_ERROR "ZeroMQ not found in ${ZEROMQ}.")
elseif(ZEROMQ_LIBRARIES)
  # ZeroMQ libraries found in system, but headers were not.
  message(STATUS "ZeroMQ headers not found. Please install the development package and the cppzmq interface. Disabling ZeroMQ support.")
else()
  message(STATUS "ZeroMQ not found. Disabling ZeroMQ support.")
endif()
