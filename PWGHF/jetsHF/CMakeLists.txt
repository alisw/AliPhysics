# Generating libPWGHFjetsHF

set(MODULE PWGHFjetsHF)

include_directories(${AliPhysics_SOURCE_DIR}/PWGHF/jetsHF)

include_directories(${ROOT_INCLUDE_DIRS}
#		    ${AliPhysics_SOURCE_DIR}/PWG/ ....
                   )

set(SRCS
   )

string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generating dictionary
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES)

# Generating the ROOT map
#  - If you are not sure how to proceed please contact aliroot-git-admins@cern.ch
set(LIBDEPS ${ROOT_DEPENDENCIES} ${ALIROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")
# End Generating the ROOT map

# Custom compilation flags
set(MODULE_COMPILE_FLAGS )
set(MODULE_LINK_FLAGS )

add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS ${MODULE_COMPILE_FLAGS})

# Additional linking flags
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS ${MODULE_LINK_FLAGS})

# Link library to dependecies
target_link_libraries(${MODULE} ${LIBDEPS})

# Library installation into $CMAKE_INSTALL_PREFIX/lib
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Header installation into $CMAKE_INSTALL_PREFIX/include
install(FILES ${HDRS} DESTINATION include)

# Macros installation
#install(DIRECTORY macros DESTINATION ${MODULE})
