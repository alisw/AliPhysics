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
    message(STATUS "Generating dictionary ${DNAME} for ${LDNAME}")
    
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
    message(STATUS "LIBNAME = ${LIBNAME}")
    message(STATUS "LIBDEPS = ${LIBDEPS}")
    message(STATUS "LINKDEF = ${LINKDEF}")
    
    message(STATUS "ROOT_LIBMAP=${ROOT_LIBMAP}")
    
    foreach(d ${LIBDEPS})
        get_filename_component(_ext ${d} EXT)
        if(_ext)
            set(Int_DEPENDENCIES ${Int_DEPENDENCIES} ${d})
        else()
            set(Int_DEPENDENCIES ${Int_DEPENDENCIES} lib${d}.so)
        endif()
    endforeach()
    
    message(STATUS "Generating ROOT map for ${LIBNAME}")
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap
                       COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_LIBMAP}
                       ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap -l lib${LIBNAME}.so -d ${Int_DEPENDENCIES} -c ${LINKDEF}
                       DEPENDS ${LIBNAME}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
                      )
    add_custom_target(lib${LIBNAME}.rootmap ALL DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap DESTINATION lib)
    
endmacro(generate_rootmap)