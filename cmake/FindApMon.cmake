# Find for ApMon library
# Setting:
#       - ApMon_LIBRARIES
#       - ApMon_INCLUDE_DIR

set(ApMon_FOUND FALSE)

if(ApMon)
    # check for the existance of the ApMon library
    find_library(ApMon_LIBRARIES NAMES apmoncpp PATHS ${ApMon}/lib NO_DEFAULT_PATH DOC "Path to ApMon library")
    
    if(NOT ApMon_LIBRARIES)
        message(FATAL_ERROR "Could not locate ApMon library inside ${ApMon}/lib")
    else()
        message(STATUS "Found ApMon library: ${ApMon_LIBRARIES}")
    endif()
    
    # check for the existance of the ApMon header
    find_path(ApMon_INCLUDE_DIR NAMES ApMon.h PATHS ${ApMon}/include NO_DEFAULT_PATH  DOC "Path to ApMon header folder.")
    
    if(NOT ApMon_INCLUDE_DIR)
        message(FATAL_ERROR "Could not find ApMon header inside ${ApMon}/include")
    else()
        message(STATUS "Found ApMon header folder: ${ApMon_INCLUDE_DIR}")
    endif()

    set(ApMon_FOUND TRUE)
    mark_as_advanced(${ApMon_LIBRARIES} ${ApMon_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Please point to the ApMon installation using -DApMon=/install/point")
endif(ApMon)