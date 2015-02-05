# Find module for DIM module
# Setting:
#       - DIM_LIBRARIES 
#       - DIM_INCLUDE_DIR

set(DIM_FOUND FASLE)

if(DIMDIR AND ODIR)
    # check for the existance of the DIM library
    find_library(DIM_LIBRARIES NAMES dim PATHS ${DIMDIR}/${ODIR}  NO_DEFAULT_PATH DOC "Path to DIM library")
    
    if(NOT DIM_LIBRARIES)
        message(FATAL_ERROR "DIM library not found inside ${DIMDIR}/${ODIR}")
    else()
        message(STATUS "Found DIM library: ${DIM_LIBRARIES}")
    endif()

    # check for the existance of the DIM header
    find_path(DIM_INCLUDE_DIR NAMES dim.h PATHS ${DIMDIR}/dim NO_DEFAULT_PATH  DOC "Path to DIM header folders.")
    
    if(NOT DIM_INCLUDE_DIR)
        message(FATAL_ERROR "DIM header not found inside ${DIM}/dim folder")
    else()
        message(STATUS "Found DIM header folder: ${DIM_INCLUDE_DIR}")
    endif()
    
    set(DIM_FOUND true)
    mark_as_advanced(DIM_LIBRARIES DIM_INCLUDE_DIR)
endif()