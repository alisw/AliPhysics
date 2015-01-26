# Checks for a ZeroMQ installation
# ZeroMQ custom installation can be pointed using -DZEROMQ
#       - ZEROMQ_INCLUDE_DIR - Where to find zeromq include sub-directory.
#       - ZEROMQ_LIBRARIES   - List of libraries when using zeromq.
#       - ZEROMQ_FOUND       - True if zeromq found.

message(STATUS "Checking for ZeroMQ ${ZEROMQ}")

set(ZEROMQ_FOUND FALSE)

if(ZEROMQ)
    # ZeroMQ is installed in a custom place
    find_library(ZEROMQ_LIBRARIES NAMES zmq
                PATHS ${ZEROMQ}/lib
                NO_DEFAULT_PATH
                DOC "Path to libzmq)"
            )

    find_path(ZEROMQ_INCLUDE_DIR NAMES zmq.h zmq_utils.h zmq.hpp}
                PATHS ${ZEROMQ}/include
                NO_DEFAULT_PATH
                DOC "Path to ZeroMQ include header files."
            )       
else(ZEROMQ)
    # Check is the library is installed on the system
    find_library(ZEROMQ_LIBRARIES NAMES zmq
                DOC "Path to libzmq)"
            )

    find_path(ZEROMQ_INCLUDE_DIR NAMES zmq.hpp zmq_utils.h}
                DOC "Path to ZeroMQ include header files."
            )
endif(ZEROMQ)

set(ZEROMQ_DISABLED FALSE)

if(NOT ZEROMQ_LIBRARIES)
    message(STATUS "ZeroMQ library not found. Disabling ZeroMQ support")
    set(ZEROMQ_DISABLED TRUE)
endif()

if(NOT ZEROMQ_INCLUDE_DIR AND NOT ZEROMQ_DISABLED)
    message(STATUS "ZeroMQ headers not found. Please install development package. Disabling ZeroMQ support")
    set(ZEROMQ_DISABLED TRUE)
endif()

if(ZEROMQ_LIBRARIES AND ZEROMQ_INCLUDE_DIR)
    message(STATUS "Found ZeroMQ ${ZEROMQ_LIBRARIES}")
    set(ZEROMQ_FOUND TRUE)
endif()




