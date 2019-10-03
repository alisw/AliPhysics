# Compiler major and minor version
#       - CLANG_MAJOR.CLANG_MINOR or
#       - GCC_MAJOR.GCC_MINOR.GCC_PATCH

message(STATUS "Found ${CMAKE_CXX_COMPILER_ID} compiler, version ${CMAKE_CXX_COMPILER_VERSION}")

# Clang compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    string(REGEX REPLACE "^.*[ ]([0-9]+)\\.[0-9].*$" "\\1" CLANG_MAJOR "${CMAKE_CXX_COMPILER_VERSION}")
    string(REGEX REPLACE "^.*[ ][0-9]+\\.([0-9]).*$" "\\1" CLANG_MINOR "${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "Compiler MAJOR ${CLANG_MAJOR}, MINOR ${CLANG_MINOR}")
endif()

# GNU compiler
if(CMAKE_COMPILER_IS_GNUCXX)
    string(REGEX REPLACE "^([0-9]+).*$"                   "\\1" GCC_MAJOR "${CMAKE_CXX_COMPILER_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*$"          "\\1" GCC_MINOR "${CMAKE_CXX_COMPILER_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" GCC_PATCH "${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "Compiler major ${GCC_MAJOR}, minor ${GCC_MINOR}, patch ${GCC_PATCH}")
endif()
