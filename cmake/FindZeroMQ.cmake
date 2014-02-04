# Copyright (c) 2010-2011 Phorm, Inc.
# License: GNU LGPL v 3.0, see http://www.gnu.org/licenses/lgpl-3.0-standalone.html 
# Author: Andrey Skryabin <andrew@zmqmessage.org>, et al.

# - Find 0mq
# Find the zeromq includes and library
#
#  ZEROMQ_INCLUDE_DIR - Where to find zeromq include sub-directory.
#  ZEROMQ_LIBRARIES   - List of libraries when using zeromq.
#  ZEROMQ_FOUND       - True if zeromq found.

IF (ZEROMQ_INCLUDE_DIR)
  # Already in cache, be silent.
  SET(ZEROMQ_FIND_QUIETLY TRUE)
ENDIF (ZEROMQ_INCLUDE_DIR)

FIND_PATH(ZEROMQ_INCLUDE_DIR zmq.hpp)

SET(ZEROMQ_NAMES zmq)
FIND_LIBRARY(ZEROMQ_LIBRARY NAMES ${ZEROMQ_NAMES} )

# Handle the QUIETLY and REQUIRED arguments and set ZEROMQ_FOUND to
# TRUE if all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  ZEROMQ DEFAULT_MSG
  ZEROMQ_LIBRARY ZEROMQ_INCLUDE_DIR 
)

IF(ZEROMQ_FOUND)
  SET( ZEROMQ_LIBRARIES ${ZEROMQ_LIBRARY} )
ELSE(ZEROMQ_FOUND)
  SET( ZEROMQ_LIBRARIES )
ENDIF(ZEROMQ_FOUND)


