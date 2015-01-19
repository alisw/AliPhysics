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
# - ROOT_INCLUDE_DIR - full path to ROOT include folder
# - ROOT_HASALIEN - ROOT was built with AliEn support
# - ROOT_HASOPENGL - ROOT was built with OpenGL support
# - ROOT_HASXML - ROOT was built with XML support
# - ROOT_FORTRAN - fortran compiler

set(ROOT_FOUND FALSE)

if(ROOTSYS)
    message(STATUS "Checking for a proper ROOT installation in ${ROOTSYS}.")

    # Setting defaults
    set(ROOT_LIBDIR ${ROOTSYS}/lib)
    set(ROOT_INCLUDE_DIR ${ROOTSYS}/include)

    # Check for root-config scripts
    find_program(ROOT_CONFIG NAMES root-config PATHS ${ROOTSYS}/bin NO_DEFAULT_PATH)

    if(NOT ROOT_CONFIG)
        message(FATAL_ERROR "Could not find root-config script.")
    endif(NOT ROOT_CONFIG)

    # Check for rlibmap
    find_program(ROOT_LIBMAP NAMES rlibmap rootcling PATHS ${ROOTSYS}/bin NO_DEFAULT_PATH)
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
    string(STRIP "${ROOT_VERSION}" ROOT_VERSION)

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
    string(STRIP "${ROOT_FEATURES}" ROOT_FEATURES)

    # Checking for ROOT libdir
    execute_process(COMMAND ${ROOT_CONFIG} --libdir OUTPUT_VARIABLE ROOT_LIBDIR ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT libdir: ${error}")
    endif(error)
    string(STRIP "${ROOT_LIBDIR}" ROOT_LIBDIR)

    # Checking for ROOT libs
    execute_process(COMMAND ${ROOT_CONFIG} --noldflags --libs OUTPUT_VARIABLE ROOT_LIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT libdir: ${error}")
    endif(error)
    string(STRIP "${ROOT_LIBS}" ROOT_LIBS)

    foreach(lib ${ROOT_LIBS})
        string(REPLACE "-rdynamic" "" new_lib ${lib})
        string(REPLACE "-l" "" lib ${new_lib})
        set(ROOT_LIBRARIES ${ROOT_LIBRARIES} ${lib})
    endforeach()
    string(STRIP "${ROOT_LIBRARIES}" ROOT_LIBRARIES)
    separate_arguments(ROOT_LIBRARIES)

    # Checking for ROOT incdir
    execute_process(COMMAND ${ROOT_CONFIG} --incdir OUTPUT_VARIABLE ROOT_INCDIR ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT incdir: ${error}")
    endif(error)
    string(STRIP "${ROOT_INCDIR}" ROOT_INCDIR)
    set(ROOT_INCLUDE_DIR ${ROOT_INCDIR})

    # Checking for glibs
    execute_process(COMMAND ${ROOT_CONFIG} --noldflags --glibs OUTPUT_VARIABLE ROOT_GLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving ROOT glibs: ${error}")
    endif(error)

    # Checking for glibs
    string(STRIP "${ROOT_GLIBS}" ROOT_GLIBS)

    foreach(lib ${ROOT_GLIBS})
        string(REPLACE "-rdynamic" "" new_lib "${lib}")
        string(REPLACE "-l" "" lib "${new_lib}")
        set(ROOT_GLIBRARIES ${ROOT_GLIBRARIES} ${lib})
    endforeach()
    string(STRIP "${ROOT_GLIBRARIES}" ROOT_GLIBRARIES)
    separate_arguments(ROOT_GLIBRARIES)

    # Checking for AliEn support
    # If AliEn support is enabled we need to point to AliEn
    execute_process(COMMAND ${ROOT_CONFIG} --has-alien OUTPUT_VARIABLE ROOT_HASALIEN ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking if ROOT was build with AliEn support: ${error}")
    endif(error)
    
    #if defined
    if(ROOT_HASALIEN)
        string(STRIP "${ROOT_HASALIEN}" ROOT_HASALIEN)
        if(ROOT_HASALIEN STREQUAL "yes")
    	    if(ALIEN)
        	add_definitions(-DWITHALIEN)
        	
        	# AliEn might bring some system libraries, we need to use them
        	if(EXISTS "${ALIEN}/lib")
        	    link_directories(${ALIEN}/lib)
        	endif()
        	
        	# api/lib should always exists
        	if(EXISTS "${ALIEN}/api/lib")
        	    link_directories(${ALIEN}/api/lib)
        	endif()
        	
        	# include for AliEn
        	if(EXISTS "${ALIEN}/include")
        	    include_directories(SYSTEM ${ALIEN}/include)
        	endif()
        	
        	# api/include always exists
        	if(EXISTS "${ALIEN}/api/include")
        	    include_directories(SYSTEM ${ALIEN}/api/include)
        	endif()
        	
        	set(ROOT_HASALIEN TRUE)
    	    else(ALIEN)
    		message(FATAL_ERROR "ROOT was build with AliEn support but no AliEn installation found. Please set \"ALIEN\" to point to your AliEn installation.")
    	    endif(ALIEN)
        else()
            set(ROOT_HASALIEN FALSE)
        endif()
    endif(ROOT_HASALIEN)

    # Checking for xml support
    execute_process(COMMAND ${ROOT_CONFIG} --has-xml OUTPUT_VARIABLE ROOT_HASXML ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking if ROOT was build with xml support: ${error}")
    endif(error)
    
    # if defined
    if(ROOT_HASXML)
        string(STRIP "${ROOT_HASXML}" ROOT_HASXML)
        if(ROOT_HASXML STREQUAL "yes")
            add_definitions(-DWITHXML)
            set(ROOT_HASXML TRUE)
        else()
            set(ROOT_HASXML FALSE)
        endif()
    endif(ROOT_HASXML)

    # Checking for OpenGL support
    execute_process(COMMAND ${ROOT_CONFIG} --has-opengl OUTPUT_VARIABLE ROOT_HASOPENGL ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking if ROOT was build with OpenGL support: ${error}")
    endif(error)
    
    if(ROOT_HASOPENGL)
        string(STRIP "${ROOT_HASOPENGL}" ROOT_HASOPENGL)
        if(ROOT_HASOPENGL STREQUAL "yes")
            set(ROOT_HASOPENGL TRUE)
        else()
            set(ROOT_HASOPENGL FALSE)
        endif()
    endif()


    # Checking for fortran compiler
    execute_process(COMMAND ${ROOT_CONFIG} --f77 OUTPUT_VARIABLE ROOT_FORTRAN ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error checking ROOT fortran compiler: ${error}")
    endif(error)
    string(STRIP "${ROOT_FORTRAN}" ROOT_FORTRAN)

    # adding the libraries and the inc dir
    link_directories(${ROOT_LIBDIR})
    include_directories(SYSTEM ${ROOT_INCLUDE_DIR})
    set(ROOT_FOUND TRUE)

    # Workaround misssing XML, VMC, Minuit from ROOT static library
    # Also adding libzma.a libfreetype.a libpcre.a
    # To be removed where https://sft.its.cern.ch/jira/browse/ROOT-6904 issue is fixed
    if(ALIROOT_STATIC)
        message(WARNING "AliRoot static build enabled! libRootExtra.a will be created.\nPlease remove when ROOT https://sft.its.cern.ch/jira/browse/ROOT-6904 issue is fixed")

        # ROOT needs xml2 but it is not build with the static version so we have to add it
        find_package(LibXml2)

        if(LIBXML2_FOUND)
        
            # adding SSL
            find_package(OpenSSL)
        
            if(OPENSSL_FOUND)
                file(GLOB _extraroot "${ROOTSYS}/montercarlo/vmc/src/*.o" "${ROOTSYS}/tree/treeplayer/src/*.o" "${ROOTSYS}/io/xmlparser/src/*.o" "${ROOTSYS}/math/minuit2/src/*.o")
                add_library(RootExtra STATIC ${_extraroot})
                set_target_properties(RootExtra PROPERTIES COMPILE_FLAGS "${LIBXML2_INCLUDE_DIR} ${OPENSSL_INCLUDE_DIR}")
                set_target_properties(RootExtra PROPERTIES LINKER_LANGUAGE CXX)
                target_link_libraries(RootExtra ${ROOT_LIBDIR}/libfreetype.a ${ROOT_LIBDIR}/libpcre.a ${ROOT_LIBDIR}/liblzma.a ${LIBXML2_LIBRARIES} ${OPENSSL_LIBRARIES})
            else()
                message(FATAL_ERROR "OpenSSL not found. Coould not generate the static RootExtra. Please install Openssl")
            endif()
        else()
            message(FATAL_ERROR "libxml2 not found. Can not generate the static RootExtra. Please install libxml2")
        endif()
    endif(ALIROOT_STATIC)
else()
    message(FATAL_ERROR "ROOT installation not found! Please point to the ROOT installation using -DROOTSYS=ROOT_INSTALL_DIR.")
endif(ROOTSYS)
