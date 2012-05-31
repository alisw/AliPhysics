# AliRoot Build System Module to find and configure ROOT
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

get_filename_component(__cmake_path ${CMAKE_COMMAND} PATH)
find_program(CPACK_COMMAND cpack ${__cmake_path})
message(STATUS "Found CPack at: ${CPACK_COMMAND}")
if(NOT CPACK_COMMAND)
  message(WARNING "CPack not found you will be unable to make packages!")
endif(NOT CPACK_COMMAND)

if(UNIX)
  set(CPACK_GENERATOR "RPM")
endif(UNIX)
#if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
#  set(CPACK_GENERATOR "TBZ2;STGZ")
#endif(CMAKE_SYSTEM_NAME STREQUAL "Linux") 
if(WIN32)
  set(CPACK_GENERATOR "NSIS;ZIP")
endif(WIN32)
add_custom_target(all_packages)

macro(add_component_package __component __targetname)
  set(PACKAGE_COMPONENT ${__component})
  set(__packageConfig CPackConfig-${PACKAGE_COMPONENT}.cmake)

  configure_file(${PROJECT_SOURCE_DIR}/cmake/CPackConfig.in ${__packageConfig} @ONLY)
  message(STATUS "Package ${__targetname} configured")
  add_custom_target(${__targetname}
    ${CPACK_COMMAND} --config "${__packageConfig}" 
#    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
  add_dependencies(all_packages ${__targetname})
endmacro(add_component_package)
