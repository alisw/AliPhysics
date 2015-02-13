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
# - DATE_LIBS - DATE libs to be linked against to reported by date-config --libs
# - DATE_LIBRARIES - DATE libs as as list as extracted from date-config --libs
# - DATE_MONLIBS - static monitorlibs reported by date-config --monitorlibs without shift
# - DATE_DYNMONLIBS - dynamic monitorlibs reported by date-config --monitorlibs without shift
# - DATE_MONLIBRARIES - DATE monitor libs as a list extracted from the date-config --monitorlibs
# - DATE_STATICMON - DATE static monitor libs needed by the DA
# - DATE_RCPROXYLIBS - rcproxylibs reported by date-config --rcproxylibs
# - DATE_RCPROXYLIBRARIES - rcproxylibs as a list as extracted from date-config --rcproxylibs

#########################
# Functions definitions
#########################

# A function to find all the libraries listed in library_list. The list contains the short
# names of libraries (ie. Db instead of libDb.so). 
# Libraries are search in library_paths. 
# It returns the list of libraries, with full path and full name.
# Author: Barthelemy Von Haller
function(find_date_libraries _output library_list library_paths)
    FOREACH (LIB ${library_list})
        # Find first static, this is used by the DA
        find_library(DATE_LIBRARY_${LIB} NAMES "lib${LIB}.a" PATHS ${library_paths})
        
        if(NOT DATE_LIBRARY_${LIB})
            find_library(DATE_LIBRARY_${LIB} NAMES ${LIB}  PATHS ${library_paths})
        endif()
        
        mark_as_advanced(DATE_LIBRARY_${LIB})
        set(_output_tmp ${_output_tmp} ${DATE_LIBRARY_${LIB}})
    ENDFOREACH (LIB ${library_list})
    
    set(${_output} ${_output_tmp} PARENT_SCOPE)
endfunction(find_date_libraries _output library_list library_paths)

# DATE_CONFIG set from the configuration
if(DATE_CONFIG)

    # Setting DIMDIR, ODIR and ROOTSYS in the environment, they are needed by date-config
    set(ENV{DIMDIR} ${DIMDIR})
    set(ENV{ODIR} ${ODIR})
    set(ENV{ROOTSYS} ${ROOTSYS})

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
    
    # DATE_LIBRARIES
    # Extracting the list of dynamic and static libraries from the --libs 
    # The list is needed during the generation of the DAs
    # using a find_library in order to get the full path and not need to use -L during linking of the DAs
    string(REGEX MATCHALL "[-]l[^- ]+" DATE_LIBRARIES_TMP ${DATE_LIBS})
    string(REGEX REPLACE "[-]l" ";" DATE_LIBRARIES_TMP ${DATE_LIBRARIES_TMP})
    # Get the list of search path using -Lyyy -> yyy
    string(REGEX MATCHALL "[-]L[^- ]+" DATE_LIBRARIES_PATH_TMP ${DATE_LIBS})
    string(REGEX REPLACE "[-]L" ";" DATE_LIBRARIES_PATH_TMP ${DATE_LIBRARIES_PATH_TMP})
    find_date_libraries(DATE_LIBRARIES "${DATE_LIBRARIES_TMP}" "${DATE_LIBRARIES_PATH_TMP}")
    
    # Fix for mysql bug https://bugs.launchpad.net/percona-server/+bug/1287374
    set(DATE_LIBS "${DATE_LIBS} -L/usr/lib64/mysql/")
    
    # setting the monlibs
    execute_process(COMMAND ${DATE_CONFIG} --monitorlibs=noshift OUTPUT_VARIABLE DATE_MONLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE monitorlibs : ${error}")
    endif(error)
    
    # If the flags are not empty we strip them
    if(DATE_MONLIBS)
        string(STRIP ${DATE_MONLIBS} DATE_MONLIBS)
    endif()

    # DATE_MONLIBRARIES
    # Extracting the list of dynamic and static libraries from the --libs 
    # The list is needed during the generation of the DAs
    # Removing everything that starts with - to leave all the static libraries
    # Replacing space with ; to create a list that will be sent to the linked 
    string(REGEX REPLACE "-[^ \r\n\t].+( |$)?" "" DATE_STATICMON ${DATE_MONLIBS})
    if(DATE_STATICMON)
        string(REPLACE " " ";"  DATE_STATICMON ${DATE_STATICMON})
    endif(DATE_STATICMON)
    
    # Extracting all the shared libraries 
    # using a find_library in order to get the full path and not need to use -L during linking of the DAs
    string(REGEX MATCHALL "[-]l[^- ]+" DATE_MONLIBRARIES_TMP ${DATE_MONLIBS})
    string(REGEX REPLACE "[-]l" ";" DATE_MONLIBRARIES_TMP ${DATE_MONLIBRARIES_TMP})
    # Get the list of search path using -Lyyy -> yyy
    string(REGEX MATCHALL "[-]L[^- ]+" DATE_MONLIBRARIES_PATH_TMP ${DATE_MONLIBS})
    string(REGEX REPLACE "[-]L" ";" DATE_MONLIBRARIES_PATH_TMP ${DATE_MONLIBRARIES_PATH_TMP})
    find_date_libraries(DATE_MONLIBRARIES "${DATE_MONLIBRARIES_TMP}" "${DATE_MONLIBRARIES_PATH_TMP}")
    set(DATE_MONLIBRARIES ${DATE_STATICMON} ${DATE_MONLIBRARIES})
    
    # Fix for mysql bug https://bugs.launchpad.net/percona-server/+bug/1287374
    set(DATE_MONLIBS "${DATE_MONLIBS} -L/usr/lib64/mysql/")

    # setting the rclibs
    execute_process(COMMAND ${DATE_CONFIG} --rcproxylibs OUTPUT_VARIABLE DATE_RCPROXYLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE rcproxylibs")
    endif(error)

    # If the flags are not empty we strip them
    if(DATE_RCPROXYLIBS)
        string(STRIP ${DATE_RCPROXYLIBS} DATE_RCPROXYLIBS)
    endif()

    # DATE_LIBRARIES
    # Extracting the list of dynamic and static libraries from the --libs 
    # The list is needed during the generation of the DAs
    string(REGEX MATCHALL "[-]l[^- ]+" DATE_RCPROXYLIBRARIES_TMP ${DATE_RCPROXYLIBS})
    string(REGEX REPLACE "[-]l" ";" DATE_RCPROXYLIBRARIES_TMP ${DATE_RCPROXYLIBRARIES_TMP})
    # Get the list of search path using -Lyyy -> yyy
    string(REGEX MATCHALL "[-]L[^- ]+" DATE_RCPROXYLIBRARIES_PATH_TMP ${DATE_RCPROXYLIBS})
    string(REGEX REPLACE "[-]L" ";" DATE_RCPROXYLIBRARIES_PATH_TMP ${DATE_RCPROXYLIBRARIES_PATH_TMP})
    find_date_libraries(DATE_RCPROXYLIBRARIES "${DATE_RCPROXYLIBRARIES_TMP}" "${DATE_RCPROXYLIBRARIES_PATH_TMP}")

    # setting the monlibs
    execute_process(COMMAND ${DATE_CONFIG} --monitorlibs=dyn OUTPUT_VARIABLE DATE_DYNMONLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving DATE monitorlibs : ${error}")
    endif(error)
    
    # If the flags are not empty we strip them
    if(DATE_DYNMONLIBS)
        string(STRIP ${DATE_DYNMONLIBS} DATE_DYNMONLIBS)
    endif()

    # unsetting all environment variables
    unset(ENV{DIMDIR})
    unset(ENV{ODIR})
    unset(ENV{ROOTSYS})

    set(DATE_FOUND TRUE)
else()
    message(STATUS "DATE not found")
endif(DATE_CONFIG)

