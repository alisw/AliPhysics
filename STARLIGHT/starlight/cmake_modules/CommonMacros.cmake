###########################################################################
#
#    Copyright 2010
#
#    This file is part of Starlight.
#
#    Starlight is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#	  
#    Starlight is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#	  
#    You should have received a copy of the GNU General Public License
#    along with Starlight. If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################
#
# File and Version Information:
# $Rev:: 27                          $: revision of last commit
# $Author:: bgrube                   $: author of last commit
# $Date:: 2010-11-29 13:00:04 +0100 #$: date of last commit
#
# Description:
#      collection of useful cmake macros
#
#
###########################################################################


# takes list of file names and returns file name list with new extension
# example:
#   switch_file_extension("${CC_LIST}" ".cc" ".h" H_LIST)
function(switch_file_extension IN_FILE_LIST OLD_EXT NEW_EXT OUT_FILE_LIST)
  set(NEW_FILE_LIST)
  foreach(_OLD_FILE ${IN_FILE_LIST})
    string(REGEX REPLACE "^(.*)${OLD_EXT}$" "\\1${NEW_EXT}" _NEW_FILE ${_OLD_FILE})
    set(NEW_FILE_LIST ${NEW_FILE_LIST} ${_NEW_FILE})
  endforeach()
  set(${OUT_FILE_LIST} ${NEW_FILE_LIST})
endfunction(switch_file_extension)


# adds standard shared library
# additional libraries that should be linked to can be given as optional arguments
function(make_shared_library LIB_NAME SOURCES)
  add_library_tested(${LIB_NAME} SHARED ${SOURCES})
  # proccess link libraries in additional arguments
  foreach(_LIB ${ARGN})
    target_link_libraries(${LIB_NAME} ${_LIB})
  endforeach()
endfunction(make_shared_library)


# adds standard executable
# additional libraries that should be linked to can be given as optional arguments
function(make_executable EXE_NAME SOURCES)
  add_executable(${EXE_NAME} ${SOURCES})
  # proccess link libraries in additional arguments
  foreach(_LIB ${ARGN})
    target_link_libraries(${EXE_NAME} ${_LIB})
  endforeach()
endfunction(make_executable)


macro(enforce_out_of_source_build)
  if(${AliRoot_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "Building this project in the source directory is not allowed. Please remove CMakeCache.txt, create a build directory, and run cmake there, for example:
rm CMakeCache.txt
mkdir build && cd build
cmake ..")
  endif()
endmacro(enforce_out_of_source_build)
