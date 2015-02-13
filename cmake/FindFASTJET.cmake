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

# AliRoot Build System Module to find and configure FASTJET
# Author: Marco van Leeuwen 
# FastJet installation is pointed during cmake configuration -DFASTJET
# Variables set during Find:
# - FASTJET_CONFIG - fastjet-config location
# - FASTJET_VERSION
# - FASTJET_VERSION_MAJOR
# - FASTJET_VERSION_MINOR
# - FASTJET_VERSION_PATCH
# - FASTJET_INCLUDE_DIR - fastjet headers location
# - FASTJET_LIBS_DIR - fastjet libraries location
# - FASTJET_DEFINITIONS - fastjet definition flags
# - FASTJET_CXXFLAGS - fastjet compilation flags
# - FASTJET_LIBS - fastjet libraries - array

set(FASTJET_FOUND FALSE)

if(FASTJET)
    # Check for fastjet-config script
    find_program(FASTJET_CONFIG NAMES fastjet-config PATHS ${FASTJET}/bin NO_DEFAULT_PATH)
    if(NOT FASTJET_CONFIG)
        message(FATAL_ERROR "Could not find fastjet-config executable")
    endif(NOT FASTJET_CONFIG)
    mark_as_advanced(FASTJET_CONFIG)

    # Check for header installation
    find_path(FASTJETHH include/fastjet/PseudoJet.hh PATHS ${FASTJET})

    if (FASTJETHH-NOTFOUND)
        message(FATAL_ERROR "Header file fastjet/PseudoJet.hh not found in ${FASTJET}/include. Please check your FASTJET installation")
    endif(FASTJETHH-NOTFOUND)
    mark_as_advanced(FASTJETHH)

    # FastJet version
    execute_process(COMMAND ${FASTJET_CONFIG} --version OUTPUT_VARIABLE FASTJET_VERSION ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving FastJet version : ${error}")
    endif(error)

    # Extract major, minor, and patch versions from
    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" FASTJET_VERSION_MAJOR "${FASTJET_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+.*" "\\1" FASTJET_VERSION_MINOR "${FASTJET_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" FASTJET_VERSION_PATCH "${FASTJET_VERSION}")
    
    if(FASTJET_VERSION_MAJOR LESS 3)
        message(FATAL_ERROR "FastJet ${FASTJET_VERSION_MAJOR}.${FASTJET_VERSION_MINOR}.${FASTJET_VERSION_PATCH} version is not suported. Please install a version >= 3.0.*")
    endif()
    
    if(FASTJET_VERSION_MINOR GREATER 0)
        message(FATAL_ERROR "FastJet ${FASTJET_VERSION_MAJOR}.${FASTJET_VERSION_MINOR}.${FASTJET_VERSION_PATCH} version is not suported.  Please install a version >= 3.0.*")
    endif()

    # Extracting compilation flags
    execute_process(COMMAND ${FASTJET_CONFIG} --cxxflags OUTPUT_VARIABLE FASTJET_CXXFLAGS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving FastJet compilation flags : ${error}")
    endif(error)
    
    # Extracting the include path from the CXX flags
    set(FASTJET_INCLUDE_DIR)
    if(FASTJET_CXXFLAGS)
	string(STRIP ${FASTJET_CXXFLAGS} FASTJET_CXXFLAGS)
	string(REGEX MATCHALL  "[-][I]([^ ])+" incfolders ${FASTJET_CXXFLAGS})
	
	foreach(incfolder ${incfolders})
	    string(REPLACE "-I" "" incfolder ${incfolder})
	    set(FASTJET_INCLUDE_DIR ${incfolder} ${FASTJET_INCLUDE_DIR})
	endforeach()
    endif(FASTJET_CXXFLAGS)

    # Extracting libraries and linking options
    execute_process(COMMAND ${FASTJET_CONFIG} --libs OUTPUT_VARIABLE FASTJET_CONFIGLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving FastJet libs : ${error}")
    endif(error)
    
    # Extracting the list of needed libraries during linking and the linking folders
    set(FASTJET_LIBS)
    set(FASTJET_LIBS_DIR)
    if(FASTJET_CONFIGLIBS)
	string(STRIP ${FASTJET_CONFIGLIBS} FASTJET_CONFIGLIBS)
    	
	# Extracting the list of libraries needed during linking
	# -l strings
	string(REGEX MATCHALL "[-][l]([^ ])+" fjlibs "${FASTJET_CONFIGLIBS}")
    
	foreach(flib ${fjlibs})
    	    string(REPLACE "-l" "" flib "${flib}")
    	    set(FASTJET_LIBS ${FASTJET_LIBS} ${flib})
	endforeach()

	# Extracting the list of libraries folders
	# -L strings
	string(REGEX MATCHALL "[-][L]([^ ])+" fjlibdirs "${FASTJET_CONFIGLIBS}")
	
	foreach(fjlibdir ${fjlibdirs})
	    string(REPLACE "-L" "" fjlibdir ${fjlibdir})
	    set(FASTJET_LIBS_DIR ${fjlibdir} ${FASTJET_LIBS_DIR})
	endforeach()
    endif()

    set(FASTJET_FOUND TRUE)
    set(FASTJET_DEFINITIONS "-DHAVE_FASTJET")
    message(STATUS "FastJet ${FASTJET_VERSION_MAJOR}.${FASTJET_VERSION_MINOR}.${FASTJET_VERSION_PATCH} installation found: ${FASTJET}")
else()
    message(STATUS "FastJet not found: disabling dependencies")
endif(FASTJET)
