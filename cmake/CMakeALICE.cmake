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

# General purpose functions

#########################
# ROOT utilities
#########################

# Generation of the dictionaries
# @DNAME  Dictionary name
# @LDNAME LinkDef file name, ex: LinkDef.h
# @DHDRS  Dictionary headers
# @DINCDIR Include folders that need to be passed to cint/cling
macro(generate_dictionary DNAME LDNAME DHDRS DINCDIRS)
    # Creating the INCLUDE path for cint/cling
    foreach( dir ${DINCDIRS})
        set(INCLUDE_PATH -I${dir} ${INCLUDE_PATH})
    endforeach()
    
    # Generate the dictionary
#    message(STATUS "Generating dictionary ${DNAME} for ${LDNAME}")
    
#    message(STATUS "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx")
#    message(STATUS "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.h")
#    message(STATUS "bbb${INCLUDE_PATH}bbb")
#    message(STATUS "${DHDRS} ${LDNAME}")
#    message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}")
    
    # Get the definitions from the directory to be sent to CINT
    get_directory_property(tmpdirdefs DEFINITIONS)
    string(REPLACE " " ";" tmpdirdefs ${tmpdirdefs})

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.h
                       COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_CINT}
                       ARGS -f ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx -c -p 
                       ${tmpdirdefs} ${INCLUDE_PATH} 
                       ${DHDRS} ${LDNAME}
                       DEPENDS ${DHDRS} ${LDNAME} ${ROOT_CINT}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      )
endmacro(generate_dictionary)

# Generate the ROOTmap files
# @LIBNAME - library name: libAnalysis.so -> Analysis.rootmap
# @LIBDEPS - library dependencies
# @LINKDEF - LinkDef header
macro(generate_rootmap LIBNAME LIBDEPS LINKDEF)
#    message(STATUS "LIBNAME = ${LIBNAME}")
#    message(STATUS "LIBDEPS = ${LIBDEPS}")
#    message(STATUS "LINKDEF = ${LINKDEF}")
#    message(STATUS "ROOT_LIBMAP=${ROOT_LIBMAP}")

    set(LOCAL_DEPS)
    foreach(file ${LIBDEPS})
        get_filename_component(ext ${file} EXT)
        if(ext)
            set(LOCAL_DEPS ${LOCAL_DEPS} ${file})
        else()
            set(LOCAL_DEPS ${LOCAL_DEPS} lib${file}.so)
        endif()
    endforeach()

#    message(STATUS "Generating ROOT map for ${LIBNAME}")
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap
                       COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_LIBMAP}
                       ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap -l lib${LIBNAME}.so -d ${LOCAL_DEPS} -c ${LINKDEF}
                       DEPENDS ${LIBNAME}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
                      )
    add_custom_target(lib${LIBNAME}.rootmap ALL DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap DESTINATION lib)
    
endmacro(generate_rootmap)

#########################
# Static utilities
#########################

# Generate the static dependecies from dynamic list
# @ shared_list - list of shared libraries
# @ static_list - the name of the variable that will contain the list of static libraries
macro(generate_static_dependencies shared_list static_list)
    set(static_list_tmp "")
    foreach(shared_lib ${shared_list})
        set(static_list_tmp ${static_list_tmp} "${shared_lib}-static")
    endforeach()
    
    # create the variable with the name received by the macro
    set(${static_list} ${static_list_tmp})
    # set the scope to parent in order to be visible in the parent
    set(${static_list} PARENT_SCOPE)
endmacro(generate_static_dependencies)

#########################
# DA utilities
#########################
# Extract DA information to be inserted into the rpm
function(getDAinfo _info _detector _daname info)
    file(STRINGS "${_detector}${_daname}da.cxx" tmpinfo REGEX "${info}:")
    string(REPLACE "${info}:\ " "" tmpinfo ${tmpinfo})
    set(${_info} ${tmpinfo} PARENT_SCOPE)
endfunction()

# DA rpm creation
macro(createDArpm DETECTOR ALGORITHM)
    getDAinfo(contact "${DETECTOR}" "${ALGORITHM}" "Contact")
    getDAinfo(link "${DETECTOR}" "${ALGORITHM}" "Link")
    getDAinfo(refrun "${DETECTOR}" "${ALGORITHM}" "Reference Run")
    getDAinfo(runtype "${DETECTOR}" "${ALGORITHM}" "Run Type")
    getDAinfo(datype "${DETECTOR}" "${ALGORITHM}" "DA Type")
    getDAinfo(evennr "${DETECTOR}" "${ALGORITHM}" "Number of events needed")
    getDAinfo(ifiles "${DETECTOR}" "${ALGORITHM}" "Input Files")
    getDAinfo(ofiles "${DETECTOR}" "${ALGORITHM}" "Output Files")
    getDAinfo(trigger "${DETECTOR}" "${ALGORITHM}" "Trigger types used")
    set(RPM_DESCRIPTION "contact: ${contact}
Link:${link}
Reference run:${refrun}
Run Type:${runtype}
DA Type:${datype}
Number of events needed: ${evennr}
Input Files:${ifiles}
Output Files:${ofiles}
Trigger types used:${trigger}")

    set(DA_EXECUTABLE "${DETECTOR}${ALGORITHM}da")
    set(DETECTOR ${DETECTOR})
    set(ALGORITHM ${ALGORITHM})
    configure_file("${AliRoot_SOURCE_DIR}/cmake/da.spec.in" "${CMAKE_CURRENT_BINARY_DIR}/${ALGORITHM}-da.spec" @ONLY)

    add_custom_command(TARGET ${DETECTOR}${ALGORITHM}da POST_BUILD
                       COMMAND mkdir ARGS -p da-${ALGORITHM}-rpm/opt/daqDA-${DETECTOR}-${ALGORITHM}/
                       COMMAND cp ARGS ${DETECTOR}${ALGORITHM}da da-${ALGORITHM}-rpm/opt/daqDA-${DETECTOR}-${ALGORITHM}/
                       COMMAND rpmbuild ARGS --verbose --define "_topdir ${CMAKE_CURRENT_BINARY_DIR}" --define "%buildroot ${CMAKE_CURRENT_BINARY_DIR}/da-${ALGORITHM}-rpm" -bb ${ALGORITHM}-da.spec
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
    )
endmacro()