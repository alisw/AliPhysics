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

###########################################################################
# ROOT utilities
###########################################################################

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

###########################################################################
# Shared librarires utilities
###########################################################################
macro(generate_shared_library)
    # Generate the dictionary
    # It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
    get_directory_property(incdirs INCLUDE_DIRECTORIES)
    set(incdirs ${MODULE_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR} ${incdirs})
    generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${MODULE_HDRS}" "${incdirs}")

    # Generate the ROOT map
    # Dependecies
    set(MODULE_LIBDEPS ${MODULE_ALIROOT_DEPENDENCIES} ${MODULE_ROOT_DEPENDENCIES})
    generate_rootmap("${MODULE}" "${MODULE_LIBDEPS}" "${MODULE}LinkDef.h")

    # Create an object to be reused in case of static libraries 
    # Otherwise the sources will be compiled twice
    add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
    # Add a library to the project using the object
    add_library(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
    target_link_libraries(${MODULE} ${MODULE_LIBDEPS})

    # Setting the correct headers for the object as gathered from the dependencies
    target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
    set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE} PUBLIC ${incdirs})

    # Setting compilation flags for the object
    set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
    # Setting the linking flags for the library
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

    # Installation
    install(TARGETS ${MODULE}
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)

    install(FILES ${MODULE_HDRS_INSTALL} DESTINATION include)
endmacro()

###########################################################################
# Static libraries utilities
###########################################################################

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

# Generate the static library
macro(generate_static_library)
    add_library(${MODULE}-static STATIC $<TARGET_OBJECTS:${MODULE}-object>)
    
    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${MODULE_ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")
    target_link_libraries(${MODULE}-static ${STATIC_ALIROOT_DEPENDENCIES} Root RootExtra)
    
    # Public include folders that will be propagated to the dependecies
    target_include_directories(${MODULE}-static PUBLIC ${incdirs})

    set_target_properties(${MODULE}-static PROPERTIES OUTPUT_NAME ${MODULE})
    set_target_properties(${MODULE}-static PROPERTIES LINK_FLAGS "-Wl,--whole-archive")

    # Installation
    install(TARGETS ${MODULE}-static
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
endmacro()

###########################################################################
# DA utilities
###########################################################################

# Extract the first comment from a DA file
# Find the position for first /* and */ and extract the substring
macro(getDAdescription _detector _daname)
    # Reading the file into a string
    file(READ "${_detector}${_daname}da.cxx" tmpinfo)
    
    # Find the first occurance of /* */
    string(FIND "${tmpinfo}" "/*" _first_position)
    string(FIND "${tmpinfo}" "*/" _second_position)
    
    # Adding and removing 2 characters to remove /* */
    math(EXPR _first_position ${_first_position}+2)
    math(EXPR _second_position ${_second_position}-2)
    
    # Generating the length of the comment in order to take out the description
    math(EXPR _desc_length ${_second_position}-${_first_position})
    
    if(${_desc_length} EQUAL 0 OR ${_desc_length} LESS 0)
        message(FATAL_ERROR "{_detector}${_daname}da.cxx does not contain a description. Please add the description as the first /*comment*/ in the file")
    else()
        string(SUBSTRING "${tmpinfo}" ${_first_position}  ${_second_position} _da_description)
        string(STRIP ${_da_description} _da_description)
        
        # The variable can be accesed by the parent
        set(RPM_DESCRIPTION ${_da_description})
    endif()
endmacro()

# Set the compilation flags
macro(setDAflags)
    # DIM
    link_directories(${DIMDIR}/${ODIR})

    #daqDA flags
    include_directories(${daqDA})
    link_directories(${daqDA})

    # AMORE definitions
    add_definitions(${AMORE_DEFINITIONS})
    include_directories(${AMORE_INCLUDE_DIR})

endmacro()

# Generate a DA
macro(generateDA DETECTOR ALGORITHM STATIC_DEPENDENCIES)
    setDAflags()

    # Generating the DA executable
    add_executable(${DETECTOR}${ALGORITHM}da ${DETECTOR}${ALGORITHM}da.cxx) #

    # DA flags and linking information
    set(MODULE_COMPILE_FLAGS)
    set(MODULE_LINK_FLAGS)

    target_link_libraries(${DETECTOR}${ALGORITHM}da ${STATIC_DEPENDENCIES} ${AMORE_AUXLIBS} daqDA ${DATE_MONLIBRARIES} ${DATE_RCPROXYLIBRARIES} Root RootExtra) # 1

    # different flags
    set(MODULE_COMPILE_FLAGS "  ${DATE_CFLAGS} ${AMORE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${AMORE_STATICLIBS}")

    set_target_properties(${DETECTOR}${ALGORITHM}da PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})
    set_target_properties(${DETECTOR}${ALGORITHM}da PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

    # Installation
    install(TARGETS ${DETECTOR}${ALGORITHM}da RUNTIME DESTINATION bin)
    
    if(DARPM)
        createDArpm("${DETECTOR}" "${ALGORITHM}")
    endif(DARPM)
endmacro()

# DA rpm creation
macro(createDArpm DETECTOR ALGORITHM)
    getDAdescription("${DETECTOR}" "${ALGORITHM}")

    set(DA_EXECUTABLE "${DETECTOR}${ALGORITHM}da")
    set(DETECTOR "${DETECTOR}")
    set(ALGORITHM "${ALGORITHM}")
    set(RPM_DESCRIPTION ${RPM_DESCRIPTION})
    
    if(ALGORITHM STREQUAL "")
        set(_ALGORITHM "none")
        set(DA_PREFIX "opt/daqDA-${DETECTOR}")
        set(DA_NAME "daqDA-${DETECTOR}")
    else()
        set(_ALGORITHM ${ALGORITHM})
        set(DA_PREFIX "opt/daqDA-${DETECTOR}-${ALGORITHM}")
        set(DA_NAME "daqDA-${DETECTOR}-${ALGORITHM}")
    endif()

    configure_file("${AliRoot_SOURCE_DIR}/cmake/da.spec.in" "${_ALGORITHM}-da.spec" @ONLY)

    add_custom_command(TARGET ${DETECTOR}${ALGORITHM}da POST_BUILD
                       COMMAND mkdir ARGS -p da-${_ALGORITHM}-rpm/root/${DA_PREFIX}/
                       COMMAND cp ARGS ${DETECTOR}${ALGORITHM}da da-${_ALGORITHM}-rpm/root/${DA_PREFIX}/
                       COMMAND rpmbuild ARGS --verbose --define "_topdir ${CMAKE_CURRENT_BINARY_DIR}/da-${_ALGORITHM}-rpm" --define "%buildroot ${CMAKE_CURRENT_BINARY_DIR}/da-${_ALGORITHM}-rpm/root" -bb ${_ALGORITHM}-da.spec
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
                       COMMENT "RPM creation for ${DETECTOR}-${_ALGORITHM}"
    )
    
    # make clean will remove also the rpm folder
    # Retrive the current list of file to be deleted - set_directory_property is overwriting, not adding to the list
    get_directory_property(_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
    set(_clean_files da-${_ALGORITHM}-rpm  ${_clean_files})
    set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${_clean_files}")
    
    # install RPM into $CMAKE_INSTALL_PREFIX/darpms
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/da-${_ALGORITHM}-rpm/RPMS/ DESTINATION darpms PATTERN "\\.rpm")
endmacro()