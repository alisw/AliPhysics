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

# Checking for a proper ROOT installation based on the ROOTSYS variable
# received during configuration
# If proper root installation it is setting the following variables
# - ROOT_VERSION - ROOT version as reported by root-config
# - ROOT_VERSION_MAJOR
# - ROOT_VERSIOM_MINOR
# - ROOT_VERSION_PATCH
# - ROOT_CONFIG - path to root-config script
# - ROOT_CINT - path to rootcint executable
# - ROOT_LIBMAP - path to rlibmap executable
# - ROOT_FEATURES - list of build features for ROOT
# - ROOT_LIBDIR - full path to ROOT library folder
# - ROOT_LIBRARIES - libraries needed for the package to be used
# - ROOT_GLIBRARIES - regular + GUI ROOT libraries + path to be used during linking
# - ROOT_INCLUDE_DIRS - full path to ROOT include folder
# - ROOT_HASALIEN - ROOT was built with AliEn support
# - ROOT_HASXML - ROOT was built with XML support
# - ROOT_FORTRAN - fortran compiler

if(ROOTSYS)
    message(STATUS "Checking for a proper ROOT installation in ${ROOTSYS}.")

    # Setting the LD_LiBRARY_PATH to point to ROOT lib folder
    set(ROOT_LIBDIR ${ROOTSYS}/lib)
    set(ROOT_INCLUDE_DIRS ${ROOTSYS}/include)

    # Check for root-config scripts
    find_program(ROOT_CONFIG NAMES root-config PATHS ${ROOTSYS}/bin NO_DEFAULT_PATH)

    if(ROOT_CONFIG)
        message(STATUS "Found root-config: ${ROOT_CONFIG}")
    else()
        message(FATAL_ERROR "Could not find root-config script.")
    endif(ROOT_CONFIG)

    # Check for rlibmap
    find_program(ROOT_LIBMAP NAMES rlibmap PATHS ${ROOTSYS}/bin NO_DEFAULT_PATH)
    if(ROOT_LIBMAP)
        message(STATUS "Found ${ROOT_LIBMAP}")
    else()
        message(FATAL_ERROR "Could not find rlibmap executable.")
    endif(ROOT_LIBMAP)

    # Check for rootcint
    find_program(ROOT_CINT NAMES rootcint PATHS ${ROOTSYS}/bin NO_DEFAULT_PATH)
    if(ROOT_CINT)
        message(STATUS "Found ${ROOT_CINT}")
    else()
        message(FATAL_ERROR "Could not find rootcint executable.")
    endif(ROOT_CINT)

    # Checking ROOT version
    execute_process(COMMAND ${ROOT_CONFIG} --version OUTPUT_VARIABLE ROOT_VERSION ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT version : ${error}")
    endif(error)
    string(STRIP ${ROOT_VERSION} ROOT_VERSION)
    
    # Extract major, minor, and patch versions from
    string(REGEX REPLACE "^([0-9]+)\\.[0-9][0-9]+\\/[0-9][0-9]+.*" "\\1" ROOT_VERSION_MAJOR "${ROOT_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9][0-9])+\\/[0-9][0-9]+.*" "\\1" ROOT_VERSION_MINOR "${ROOT_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9][0-9]+\\/([0-9][0-9]+).*" "\\1" ROOT_VERSION_PATCH "${ROOT_VERSION}")
    message(STATUS "Found ROOT version ${ROOT_VERSION_MAJOR}.${ROOT_VERSION_MINOR}.${ROOT_VERSION_PATCH}")

    # Print ROOT features
    execute_process(COMMAND ${ROOT_CONFIG} --features OUTPUT_VARIABLE ROOT_FEATURES ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT features : ${error}")
    else()
        message(STATUS "ROOT was build with the following features: ${ROOT_FEATURES}")
    endif(error)
    string(STRIP ${ROOT_FEATURES} ROOT_FEATURES)
    
    # Checking for ROOT libdir
    execute_process(COMMAND ${ROOT_CONFIG} --libdir OUTPUT_VARIABLE ROOT_LIBDIR ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT libdir: ${error}")
    endif(error)
    string(STRIP ${ROOT_LIBDIR} ROOT_LIBDIR)

    # Checking for ROOT libs
    execute_process(COMMAND ${ROOT_CONFIG} --noldflags --libs OUTPUT_VARIABLE ROOT_LIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT libdir: ${error}")
    endif(error)
    string(STRIP ${ROOT_LIBS} ROOT_LIBS)
    
    foreach(lib ${ROOT_LIBS})
        string(REPLACE "-rdynamic" "" new_lib ${lib})
        string(REPLACE "-l" "" lib ${new_lib})
        set(ROOT_LIBRARIES ${ROOT_LIBRARIES} ${lib})
    endforeach()
    string(STRIP ${ROOT_LIBRARIES} ROOT_LIBRARIES)
    separate_arguments(ROOT_LIBRARIES)
    
    # Checking for ROOT incdir
    execute_process(COMMAND ${ROOT_CONFIG} --incdir OUTPUT_VARIABLE ROOT_INCDIR ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT incdir: ${error}")
    endif(error)
    string(STRIP ${ROOT_INCDIR} ROOT_INCDIR)
    
    # Checking for glibs
    execute_process(COMMAND ${ROOT_CONFIG} --noldflags --glibs OUTPUT_VARIABLE ROOT_GLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT glibs: ${error}")
    endif(error)

    # Checking for glibs
    string(STRIP ${ROOT_GLIBS} ROOT_GLIBS)

    foreach(lib ${ROOT_GLIBS})
        string(REPLACE "-rdynamic" "" new_lib "${lib}")
        string(REPLACE "-l" "" lib "${new_lib}")
        set(ROOT_GLIBRARIES ${ROOT_GLIBRARIES} ${lib})
    endforeach()
    string(STRIP ${ROOT_GLIBRARIES} ROOT_GLIBRARIES)
    separate_arguments(ROOT_GLIBRARIES)
    
    # Checking for AliEn support
    execute_process(COMMAND ${ROOT_CONFIG} --has-alien OUTPUT_VARIABLE ROOT_HASALIEN ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking if ROOT was build with AliEn support: ${error}")
    endif(error)
    string(STRIP ${ROOT_HASALIEN} ROOT_HASALIEN)
    
    # Checking for xml support
    execute_process(COMMAND ${ROOT_CONFIG} --has-xml OUTPUT_VARIABLE ROOT_HASXML ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking if ROOT was build with xml support: ${error}")
    endif(error)
    string(STRIP ${ROOT_HASXML} ROOT_HASXML)


    # Checking for xml support
    execute_process(COMMAND ${ROOT_CONFIG} --f77 OUTPUT_VARIABLE ROOT_FORTRAN ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking ROOT fortran compiler: ${error}")
    endif(error)
    string(STRIP ${ROOT_FORTRAN} ROOT_FORTRAN)
    
    set(ROOT_FOUND TRUE)
else()
    message(FATAL_ERROR "ROOT installation not found! Please point to the ROOT installation using -DROOTSYS=ROOT_INSTALL_DIR.")
endif(ROOTSYS)