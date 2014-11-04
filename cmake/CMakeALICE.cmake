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
    
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.h
                       COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_CINT}
                       ARGS -f ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx -c -p 
                       ${INCLUDE_PATH} 
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