# -*- mode: cmake -*-

#--------------------------------------------------------------------------------#
# Linuxx8664gcc CMake System configuration file for the AliRoot Build System     #
#                                                                                # 
# Author: Johny Jose (johny.jose@cern.ch)                                        #
#                                                                                #
#--------------------------------------------------------------------------------#


cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#Retrieve Compiler Version


execute_process (COMMAND ${CMAKE_C_COMPILER} -dumpversion | cut -d. -f1 
                 OUTPUT_VARIABLE CCMAJORV
                 OUTPUT_STRIP_TRAILING_WHITESPACE) 
execute_process (COMMAND ${CMAKE_C_COMPILER} -dumpversion | cut -d. -f2
                 OUTPUT_VARIABLE CCMINORV
                 OUTPUT_STRIP_TRAILING_WHITESPACE)


# Global Optimization
set(OPT)

# ------- Setting optimization flags for default configuration -------

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))
    set(DEFAULT_CXX_FLAGS "-O -g")
    set(OPT "${DEFAULT_CXX_FLAGS}")
    message("-- Setting compiler flags for default configuration: ${DEFAULT_CXX_FLAGS}")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))

# --------------------------------------------------------------------


set(NOOPT "-g")

set(CXXOPT ${OPT})
set(CXXNOPT ${NOOPT})
set(COPT ${OPT})
set(FOPT ${OPT})

#CERNLIB defines
set(CLIBDEFS "-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ")
set(CLIBCXXOPTS)
set(CLIBCOPT)
set(CLIBFOPT ${CLIBDEFS})

set(CXXWARN "-Wall -Wno-long-long -W -Weffc++ -Wshadow -Woverloaded-virtual -ansi")

if(CCMAJORV STREQUAL "2")

  set(CXXFLAGS "${OPT} -fPIC -pipe")
  set(CXXFLAGSNO "${NOOPT} -fPIC -pipe")

elseif(CCMAJORV STREQUAL "3")

  set(CXXFLAGS "${OPT} -fPIC -pipe -fmessage-length=0 -Dlinux")
  add_definitions(-Dlinux)
  set(CXXFLAGSNO "${NOOPT} -fPIC -pipe -fmessage-length=0")

elseif(CCMAJORV STREQUAL "4")
  
  set(CXXFLAGS "${OPT} -fPIC -pipe -fmessage-length=0 -Dlinux")
  add_definitions(-Dlinux)
  set(CXXFLAGSNO "${NOOPT} -fPIC -pipe -fmessage-length=0")

else ()

  set(CXXFLAGS "${OPT} -pipe -fmessage-length=0 -Dlinux")
  add_definitions(-Dlinux)
  set(CXXFLAGSNO "${NOOPT} -pipe -fmessage-length=0") 

endif(CCMAJORV STREQUAL "2")

set(CFLAGS	"${OPT} -Wall -Werror -fPIC -pipe -Wno-long-long -pedantic-errors -ansi")
set(FFLAGS "${CLIBFOPT} ${FOPT} -fPIC -fno-second-underscore -fPIC -fno-f2c")

execute_process(COMMAND find /lib64 -name 'libNoVersion*.so' 
                COMMAND xargs --replace basename {} .so 
                COMMAND sed -e "s/lib64/ -l/"
                OUTPUT_VARIABLE LIBNOVER
                OUTPUT_STRIP_TRAILING_WHITESPACE)
                

#set(SYSLIBS "-ldl -lcrypt -L/usr/X11R6/lib -lX11 -lGL -lGLU ${LIBNOVER}")
set(SYSLIBS "-ldl -lcrypt -L/usr/X11R6/lib -lX11 ${LIBNOVER}")

if(${CMAKE_Fortran_COMPILER} MATCHES "g95")
  
  add_definitions(-DFORTRAN_G95)
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --print-search-dirs 
                  OUTPUT_VARIABLE SHLIB
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX MATCH "^.*install:[^\n]*" SHLIB ${SHLIB})
  string(REGEX REPLACE "^.*install: " "" SHLIB ${SHLIB})
  set(SHLIB "-L${SHLIB} -lf95")

elseif(${CMAKE_Fortran_COMPILER} MATCHES "gfortran")
  
  set(FFLAGS "-DFORTRAN_GFORTRAN ${FFLAGS}")
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.so
                  OUTPUT_VARIABLE _shlib
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortranbegin.a
                  OUTPUT_VARIABLE SHLIB
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(SHLIB "${_shlib} ${SHLIB}")
  set(SYSLIBS "${SYSLIBS} ${SHLIB}")

else()

  set(SHLIB "-lg2c")
  set(SYSLIBS "-lg2c")
  
endif(${CMAKE_Fortran_COMPILER} MATCHES "g95")

set(LDFLAGS "${OPT}")
set(SOFLAGS "${OPT} -shared -Wl")
set(ALLIB)



