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
# Authors: Marco van Leeuwen, Dario Berzano
# FastJet installation is pointed during cmake configuration -DFASTJET
#
# Variables set during Find:
#
# - FASTJET_CONFIG - fastjet-config location
# - FASTJET_VERSION
# - FASTJET_VERSION_MAJOR
# - FASTJET_VERSION_MINOR
# - FASTJET_VERSION_PATCH
# - FASTJET_INCLUDE_DIR - fastjet headers location
# - FASTJET_LIBS_DIR - fastjet libraries location
# - FASTJET_DEFINITIONS - fastjet definition flags
# - FASTJET_CXXFLAGS - fastjet compilation flags
# - FASTJET_LIBS - fastjet libraries (list)
# - FASTJET_ROOTDICT_OPTS - fastjet options to generate ROOT 6+ dictionaries

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
        message(FATAL_ERROR "Header file fastjet/PseudoJet.hh not found in ${FASTJET}/include. Please check your FastJet installation")
    endif(FASTJETHH-NOTFOUND)
    mark_as_advanced(FASTJETHH)

    # FastJet version
    execute_process(COMMAND ${FASTJET_CONFIG} --version OUTPUT_VARIABLE FASTJET_VERSION ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(error)
        message(FATAL_ERROR "Error retrieving FastJet version: ${error}")
    endif(error)

    # Extract major, minor, and patch versions
    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" FASTJET_VERSION_MAJOR "${FASTJET_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+.*" "\\1" FASTJET_VERSION_MINOR "${FASTJET_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" FASTJET_VERSION_PATCH "${FASTJET_VERSION}")

    if(FASTJET_VERSION_MAJOR LESS 3)
        message(FATAL_ERROR "FastJet ${FASTJET_VERSION_MAJOR}.${FASTJET_VERSION_MINOR}.${FASTJET_VERSION_PATCH} not suported: install at least >= 3.0.*")
    endif()

    # Extracting compilation flags
    execute_process(COMMAND ${FASTJET_CONFIG} --cxxflags OUTPUT_VARIABLE FASTJET_CXXFLAGS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(error)
        message(FATAL_ERROR "Error retrieving FastJet compilation flags: ${error}")
    endif(error)

    # Extracting the include path(s?) from the CXX flags (LIST)
    set(FASTJET_INCLUDE_DIR)
    if(FASTJET_CXXFLAGS)
        string(REGEX MATCHALL "(^| )-I[^ ]+" incfolders ${FASTJET_CXXFLAGS})
        foreach(incfolder ${incfolders})
            string(STRIP ${incfolder} incfolder)
            string(SUBSTRING ${incfolder} 2 -1 incfolder)
            list(APPEND FASTJET_INCLUDE_DIR ${incfolder})
        endforeach()
    endif(FASTJET_CXXFLAGS)

    # Extracting libraries and linking options
    execute_process(COMMAND ${FASTJET_CONFIG} --libs
                    OUTPUT_VARIABLE FASTJET_CONFIGLIBS
                    ERROR_VARIABLE error
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(error)
        message(FATAL_ERROR "Error retrieving FastJet libs: ${error}")
    endif(error)

    # Extracting the list of needed libraries during linking and the linking directory
    set(FASTJET_LIBS)
    set(FASTJET_LIBS_DIR)
    if(FASTJET_CONFIGLIBS)

        # Extract libraries needed at link time (-l) (LIST)
        string(REGEX MATCHALL "(^| )-l[^ ]+" fjlibs ${FASTJET_CONFIGLIBS})
        foreach(fjlib ${fjlibs})
            string(STRIP ${fjlib} fjlib)
            string(SUBSTRING ${fjlib} 2 -1 fjlib)
            list(APPEND FASTJET_LIBS ${fjlib})
        endforeach()

        # CGAL for FastJet
        list(FIND FASTJET_LIBS "CGAL" fastjetcgal)
        if(fastjetcgal GREATER -1)
            message(STATUS "FastJet was compiled with CGAL")
            # FastJet was compiled with CGAL. We need CGAL.
            if(CGAL)
                find_library(CGAL_LIB "CGAL" "${CGAL}/lib" NO_DEFAULT_PATH)
                if(NOT CGAL_LIB)
                    message(FATAL_ERROR "CGAL not found under the given location")
                endif()
                message(STATUS "CGAL found in ${CGAL}")
                list(APPEND FASTJET_LIBS_DIR ${CGAL}/lib)
            else(CGAL)
                # Attempt to use system's CGAL
                find_library(CGAL_LIB "CGAL")
                if(NOT CGAL_LIB)
                    message(FATAL_ERROR "CGAL not installed: install it or specify a custom one with -DCGAL")
                endif()
            endif(CGAL)
            list(APPEND fjextralibs "CGAL")
        else()
            message(STATUS "FastJet was compiled witouth CGAL")
        endif(fastjetcgal GREATER -1)

        # GMP for FastJet
        list(FIND FASTJET_LIBS "gmp" fastjetgmp)
        if(fastjetgmp GREATER -1)
            message(STATUS "FastJet was compiled with GMP")
            # FastJet was compiled with GMP. We need GMP.
            if(GMP)
                find_library(GMP_LIB "gmpxx" "${GMP}/lib" NO_DEFAULT_PATH)
                if(GMP_LIB MATCHES "\\.a$")
                    message(STATUS "GMP (static) found in ${GMP}: not adding flags")
                    list(REMOVE_ITEM FASTJET_LIBS "gmp")
                elseif(NOT GMP_LIB)
                    message(FATAL_ERROR "GMP not found under the given location")
                else()
                    message(STATUS "GMP found in ${GMP}")
                    list(APPEND fjextralibs "gmp")
                    list(APPEND FASTJET_LIBS_DIR ${GMP}/lib)
                endif()
            else(GMP)
                # Attempt to use system's GMP
                find_library(GMP_LIB "gmpxx")
                if(NOT GMP_LIB)
                    message(FATAL_ERROR "GMP not installed: install it or specify a custom one with -DGMP")
                endif()
                list(APPEND fjextralibs "gmp")
            endif(GMP)
        else()
            message(STATUS "FastJet was compiled witouth GMP")
        endif(fastjetgmp GREATER -1)

        # MPFR for FastJet
        list(FIND FASTJET_LIBS "mpfr" fastjetmpfr)
        if(fastjetmpfr GREATER -1)
            message(STATUS "FastJet was compiled with MPFR")
            # FastJet was compiled with MPFR. We need MPFR.
            if(MPFR)
                find_library(MPFR_LIB "mpfr" "${MPFR}/lib" NO_DEFAULT_PATH)
                if(MPFR_LIB MATCHES "\\.a$")
                    message(STATUS "MPFR (static) found in ${MPFR}: not adding flags")
                    list(REMOVE_ITEM FASTJET_LIBS "mpfr")
                elseif(NOT MPFR_LIB)
                    message(FATAL_ERROR "MPFR not found under the given location")
                else()
                    message(STATUS "MPFR found in ${MPFR}")
                    list(APPEND fjextralibs "mpfr")
                    list(APPEND FASTJET_LIBS_DIR ${MPFR}/lib)
                endif()
            else(MPFR)
                # Attempt to use system's MPFR
                find_library(MPFR_LIB "mpfr")
                if(NOT MPFR_LIB)
                    message(FATAL_ERROR "MPFR not installed: install it or specify a custom one with -DMPFR")
                endif()
                list(APPEND fjextralibs "mpfr")
            endif(MPFR)
        else()
            message(STATUS "FastJet was compiled witouth MPFR")
        endif(fastjetmpfr GREATER -1)

        # Extract library paths used by FastJet to find needed libs: (LIST)
        string(REGEX MATCHALL "(^| )-L[^ ]+" fjlibdirs ${FASTJET_CONFIGLIBS})
        foreach(fjlibdir ${fjlibdirs})
            string(STRIP ${fjlibdir} fjlibdir)
            string(SUBSTRING ${fjlibdir} 2 -1 fjlibdir)
            list(APPEND FASTJET_LIBS_DIR ${fjlibdir})
        endforeach()

        # Check if required FJ and FJContrib libs can be found
        set(fjlibs "fastjetcontribfragile;fastjetplugins;siscone_spherical;siscone;fastjettools;fastjet")
        foreach(fjlib ${fjlibs})
            set(fjl fjl-NOTFOUND)
            find_library(fjl NAMES ${fjlib}
                             PATHS ${FASTJET_LIBS_DIR} NO_DEFAULT_PATH)
            if(NOT fjl)
                message(FATAL_ERROR "Required FastJet library ${fjlib} not found!")
            endif()
        endforeach()

    endif(FASTJET_CONFIGLIBS)

    set(FASTJET_FOUND TRUE)
    set(FASTJET_DEFINITIONS "-DHAVE_FASTJET")

    if(ROOT_VERSION_MAJOR GREATER 5)
        # ROOT 6+: additional options required for dictionary generation.
        set(FASTJET_ROOTDICT_OPTS "-I${FASTJET_INCLUDE_DIR}")
        foreach(fjl ${fjlibs};${fjextralibs})
            set(FASTJET_ROOTDICT_OPTS "${FASTJET_ROOTDICT_OPTS} -rml lib${fjl}")
        endforeach()
    else(ROOT_VERSION_MAJOR GREATER 5)
        if(NOT FASTJET_VERSION VERSION_LESS "3.2.0")
            # Disable some C++11 constructs on ROOT 5 and FastJet >= 3.2.0
            set(FASTJET_ROOTDICT_OPTS "${FASTJET_ROOTDICT_OPTS} -Doverride=")
        endif()
    endif(ROOT_VERSION_MAJOR GREATER 5)

    message(STATUS "FastJet options for ROOT dictionary generation: ${FASTJET_ROOTDICT_OPTS}")
    message(STATUS "FastJet ${FASTJET_VERSION_MAJOR}.${FASTJET_VERSION_MINOR}.${FASTJET_VERSION_PATCH} installation found: ${FASTJET}")
else()
    message(STATUS "FastJet not found: disabling dependencies")
endif(FASTJET)
