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

# Find DATE DAQ installation using date-config
# Requires: 
#       - date-config has to be found in the PATH
#       - Environment variables
#               - DATE_ROOT // try to not use
#               - DATE_COMMON_DEFS // try to not use
#               - DATE_MONITOR_DIR //try not to use
# - DATE_VERSION - DATE version as reported by date-config
# - DATE_VERSION_MAJOR
# - DATE_VERSIOM_MINOR
# - DATE_CONFIG - path to date-config script
# - DATE_CFLAGS - cflags reported by date-config
# - DATE_LDFLAGS - ldflags reported by date-config
# - DATE_LIBS - DATE libs to be linked against to
# - DATE_MONLIBS - monitorlibs reported by date-config

find_program(DATE_CONFIG date-config)

if(DATE_CONFIG)
    # Checking DATE version
    execute_process(COMMAND ${DATE_CONFIG} --version OUTPUT_VARIABLE DATE_VERSION ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE version : ${error}")
    endif(error)
    string(STRIP ${DATE_VERSION} DATE_VERSION)

    # Extract major, minor, and patch versions from
    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+" "\\1" DATE_VERSION_MAJOR "${DATE_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9]+)" "\\1" DATE_VERSION_MINOR "${DATE_VERSION}")
    message(STATUS "DATE version ${DATE_VERSION_MAJOR}.${DATE_VERSION_MINOR} found.")
    
    # Checking if the environment is properly set
#    if(NOT DEFINED ENV{DATE_ROOT})
#        message(FATAL_ERROR "date-config found. Please set DATE_ROOT environment variable")
#    else()
#        set(DATE_ROOT ENV{DATE_ROOT})
#        message(STATUS "DATE_ROOT found ${DATE_ROOT}")
#    endif()

#    if(NOT DEFINED ENV{DATE_COMMON_DEFS})
#        message(FATAL_ERROR "date-config found. Please set DATE_COMMON_DEFS environment variable")
#    else()
#        set(DATE_COMMON_DEFS ENV{DATE_COMMON_DEFS})
#        message(STATUS "DATE_COMMON_DEFS found ${DATE_COMMON_DEFS}")
#    endif()

#    if(NOT DEFINED ENV{DATE_MONITOR_DIR})
#        message(FATAL_ERROR "date-config found. Please set DATE_MONITOR_DIR environment variable")
#    else()
#        set(DATE_MONITOR_DIR ENV{DATE_MONITOR_DIR})
#        message(STATUS "DATE_MONITOR_DIR found ${DATE_MONITOR_DIR}")
#    endif()

    # setting the cflags
    execute_process(COMMAND ${DATE_CONFIG} --cflags OUTPUT_VARIABLE DATE_CFLAGS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE cflags : ${error}")
    endif(error)
    
    # If flags not empty we strip them
    if(DATE_CFLAGS)
        string(STRIP ${DATE_CFLAGS} DATE_CFLAGS)
    endif()
    
    set(DATE_CFLAGS "-DALI_DATE ${DATE_CFLAGS}")

    # setting the ldflags
    execute_process(COMMAND ${DATE_CONFIG} --ldflags OUTPUT_VARIABLE DATE_LDFLAGS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE ldflags : ${error}")
    endif(error)
    
    # If the flags are not empty we strip them
    if(DATE_LDFLAGS)
        string(STRIP ${DATE_LDFLAGS} DATE_LDFLAGS)
    endif()

    # setting the libs
    execute_process(COMMAND ${DATE_CONFIG} --libs OUTPUT_VARIABLE DATE_LIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE libs : ${error}")
    endif(error)
    
    # If the flags are not empty we strip them
    if(DATE_LIBS)
        string(STRIP ${DATE_LIBS} DATE_LIBS)
    endif()
    
    # setting the monlibs
    execute_process(COMMAND ${DATE_CONFIG} --monitorlibs OUTPUT_VARIABLE DATE_MONLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE monitorlibs : ${error}")
    endif(error)
    
    # If the flags are not empty we strip them
    if(DATE_MONLIBS)
        string(STRIP ${DATE_MONLIBS} DATE_MONLIBS)
    endif()

    set(DATE_FOUND TRUE)
else()
    message(STATUS "DATE not found")
    
#    set(DATEFLAGS "-D${CMAKE_SYSTEM_NAME} -DDATE_SYS=${CMAKE_SYSTEM_NAME} -Dlong32='int' -Dlong64='long long' -DdatePointer='long'")
endif(DATE_CONFIG)

