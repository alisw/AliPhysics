# General purpose functions

#########################
# ROOT utilities
#########################

# Generation of the dictionaries
# @DNAME  Dictionary name
# @LDNAME LinkDef file name, ex: LinkDef.h
# @DHDRS  Dictionary headers
# @DINCDIR Include folders that need to be passed to cint/cling
# @EXTRADEFINITIONS - optional, extra compile flags specific to library
#       - used as ${ARGV4}
macro(generate_dictionary DNAME LDNAME DHDRS DINCDIRS)

    # Creating the INCLUDE path for cint/cling
    foreach(dir ${DINCDIRS})
        set(INCLUDE_PATH -I${dir} ${INCLUDE_PATH})
    endforeach()

    # Get the list of definitions from the directory to be sent to CINT
    get_directory_property(tmpdirdefs COMPILE_DEFINITIONS)
    foreach(dirdef ${tmpdirdefs})
        string(REPLACE "\"" "\\\"" dirdef_esc ${dirdef})
        set(GLOBALDEFINITIONS -D${dirdef_esc} ${GLOBALDEFINITIONS})
    endforeach()

    # Custom definitions specific to library
    # Received as the forth optional argument
    separate_arguments(EXTRADEFINITIONS UNIX_COMMAND "${ARGV4}")

    if(ROOT_VERSION_MAJOR LESS 6)
        # ROOT version < 6
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.h
                           COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_CINT}
                           ARGS -f ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx -c -p
                           ${GLOBALDEFINITIONS} ${EXTRADEFINITIONS} ${INCLUDE_PATH}
                           ${DHDRS} ${LDNAME}
                           DEPENDS ${DHDRS} ${LDNAME} ${ROOT_CINT}
                           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else()
      # ROOT version >= 6
      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}_rdict.pcm
                         COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_CINT}
                         ARGS -f ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx
                              -rmf ${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap
                              -rml lib${DNAME}
                              ${GLOBALDEFINITIONS}
                              ${EXTRADEFINITIONS}
                              ${INCLUDE_PATH}
                              ${DHDRS}
                              ${LDNAME}
                         DEPENDS ${DHDRS} ${LDNAME} ${ROOT_CINT}
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap" DESTINATION lib)
      install(FILES "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}_rdict.pcm" DESTINATION lib)
    endif()

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

if (ROOT_VERSION_MAJOR LESS 6)

    set(LOCAL_DEPS)
    foreach(file ${LIBDEPS})
        get_filename_component(ext ${file} EXT)
        if(ext)
            set(LOCAL_DEPS ${LOCAL_DEPS} ${file})
        else()
            set(LOCAL_DEPS ${LOCAL_DEPS} lib${file})
        endif()
    endforeach()

#    message(STATUS "Generating ROOT map for ${LIBNAME}")
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap
                       COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_LIBMAP}
                       ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap -l lib${LIBNAME} -d ${LOCAL_DEPS} -c ${LINKDEF}
                       DEPENDS ${LIBNAME}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
                      )
    add_custom_target(lib${LIBNAME}.rootmap ALL DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap DESTINATION lib)

endif (ROOT_VERSION_MAJOR LESS 6)

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


# Prepend prefix to every element in the list. Note: this function modifies the input variable: this
# does not work for macros in CMake, only for functions. Also note that it does NOT automatically
# add a / between prefix and list item as it does not assume that we are dealing with directories
function(prepend_prefix INLIST PREFIX)
    foreach(_ITEM ${${INLIST}})
        list(APPEND _OUTLIST ${PREFIX}${_ITEM})
    endforeach()
    set(${INLIST} ${_OUTLIST} PARENT_SCOPE)
endfunction()


# This function is a drop-in replacement for the following CMake command:
#
#   install(FILES ... DESTINATION ... [OTHER_ARGS])
#
# The above command takes every single file and puts it in the destination directory, but relative
# paths are not taken into consideration, i.e. files a/b/c/file.h and boo.h will be both installed
# in dest.
#
# By replacing install() with install_relative(), boo.h will end in dest, and a/b/c/file.h will end
# in dest/a/b/c/file.h, i.e. relative paths are taken into account.
#
# If an absolute path was specified for an input file, a fatal error will be raised: only relative
# paths can be specified.
#
# Since it is a drop-in command, its syntax is identical to install():
#
#   install_relative(FILES ... DESTINATION ... [OTHER_ARGS])
#
# where OTHER_ARGS is passed as-is to the underlying install() command
function(install_relative)

    set(_EXPECT_FILE TRUE)
    set(_EXPECT_DEST FALSE)
    set(_EXPECT_REST FALSE)

    foreach(_ARG ${ARGN})

        if(_EXPECT_FILE)

            if(${_ARG} STREQUAL "FILES")
                set(_EXPECT_FILE FALSE)
            else()
                message(FATAL_ERROR "You may only use install_relative() in place of install(FILES ...)")
            endif()

        elseif(_EXPECT_REST)
            # Remaining arguments
            list(APPEND _REST ${_ARG})
        elseif(_EXPECT_DEST)
            # Destination prefix
            set(_DEST ${_ARG})
            set(_EXPECT_DEST FALSE)
            set(_EXPECT_REST TRUE)
        elseif(_ARG STREQUAL "DESTINATION")
            # From now on, copy the arguments ditto to the install() command
            set(_EXPECT_DEST TRUE)
        else()
            # Append files to install
            list(APPEND _FILES ${_ARG})
        endif()

    endforeach()

    # Print out our results (debug)
    #message(STATUS "[install_relative] FILES: ${_FILES}")
    #message(STATUS "[install_relative] DEST: ${_DEST}")
    #message(STATUS "[install_relative] REST: ${_REST}")

    # Prepare a distinct install command for each file, depending on its path
    foreach(_FILE ${_FILES})

        if(CMAKE_VERSION VERSION_LESS "2.8.12")
            get_filename_component(_FILEPREFIX ${_FILE} PATH)
        else()
            get_filename_component(_FILEPREFIX ${_FILE} DIRECTORY)
        endif()
        #message(STATUS "[install_relative] ${_FILE} --> ${_FILEPREFIX}")

        string(SUBSTRING ${_FILE} 0 1 _FILE_FIRST)
        if(${_FILE_FIRST} STREQUAL "/")
            # An absolute path was found: not supported, error
            message(FATAL_ERROR "Absolute paths are not supported by install_relative(): ${_FILE}")
        endif()

        install(FILES ${_FILE} DESTINATION ${_DEST}/${_FILEPREFIX} ${_REST})

    endforeach()

endfunction()
