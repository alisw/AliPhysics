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

# Module
set(MODULE alimdc)
add_definitions(-D_MODULE_="${MODULE}")
include_directories(${AliRoot_SOURCE_DIR}/${MODULE})

# Aditional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR} 
                    ${AliRoot_SOURCE_DIR}/RAW/MDC
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT alimdc_main.cxx)
# Add a library to the project using the object
add_executable(${MODULE} $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} STEERBase EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Physics VMC Thread ESD MDC RAWDatabase RAWDatarec)
 
# Compile and link flags
set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
endif(DATE_FOUND)

# Setting target compile and linking flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE} RUNTIME DESTINATION bin)

# Creating the AliMDC static library and the rpm
if(MDCRPM)
    add_library(AliMDC STATIC $<TARGET_OBJECTS:${MODULE}-object> $<TARGET_OBJECTS:MDC-object> $<TARGET_OBJECTS:RAWDatarec-object> $<TARGET_OBJECTS:RAWDatabase-object> $<TARGET_OBJECTS:ESD-object> $<TARGET_OBJECTS:STEERBase-object>)

    # list of shared dependencies / the name of the variable containing the list of static ones
    generate_static_dependencies("${ALIROOT_DEPENDENCIES}" "STATIC_ALIROOT_DEPENDENCIES")
    target_link_libraries(AliMDC ${STATIC_ALIROOT_DEPENDENCIES} Root RootExtra)

    set_target_properties(AliMDC PROPERTIES LINK_FLAGS "-Wl,--whole-archive -Wl,--no-as-needed")

    # Installation
    install(TARGETS AliMDC
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)

    configure_file("alimdc.spec.in" "alimdc.spec" @ONLY)

    add_custom_command(TARGET alimdc POST_BUILD
                       COMMAND mkdir ARGS -p alimdc-rpm/root/opt/alimdc/lib
                       COMMAND mkdir ARGS -p alimdc-rpm/root/opt/alimdc/include
                       COMMAND cp ARGS libAliMDC.a alimdc-rpm/root/opt/alimdc/lib
                       COMMAND cp ARGS ${ROOT_LIBDIR}/libRoot.a alimdc-rpm/root/opt/alimdc/lib/
                       COMMAND cp ARGS ${CMAKE_BINARY_DIR}/libRootExtra.a alimdc-rpm/root/opt/alimdc/lib/
                       COMMAND cp ARGS ${ROOT_LIBDIR}/libpcre.a alimdc-rpm/root/opt/alimdc/lib/
                       COMMAND cp ARGS ${ROOT_LIBDIR}/liblzma.a alimdc-rpm/root/opt/alimdc/lib/
                       COMMAND cp ARGS ${ROOT_LIBDIR}/libfreetype.a alimdc-rpm/root/opt/alimdc/lib/
                       COMMAND cp ARGS ${CMAKE_BINARY_DIR}/libRootExtra.a alimdc-rpm/root/opt/alimdc/lib/
                       COMMAND cp ARGS ${AliRoot_SOURCE_DIR}/RAW/MDC/mdc.h alimdc-rpm/root/opt/alimdc/include
                       COMMAND rpmbuild ARGS --verbose --define "_topdir ${CMAKE_CURRENT_BINARY_DIR}/alimdc-rpm" --define "%buildroot ${CMAKE_CURRENT_BINARY_DIR}/alimdc-rpm/root" -bb alimdc.spec
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
                       COMMENT "RPM creation for AliMDC"
    )
    
    # make clean will remove also the rpm folder
    # Retrive the current list of file to be deleted - set_directory_property is overwriting, not adding to the list
    get_directory_property(_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
    set(_clean_files alimdc-rpm ${_clean_files})
    set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${_clean_files}")
    
    # install RPM into $CMAKE_INSTALL_PREFIX/darpms
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/alimdc-rpm/RPMS/ DESTINATION mdcrpm PATTERN "\\.rpm")

endif(MDCRPM)
