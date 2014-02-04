# -*- mode: cmake -*-

#--------------------------------------------------------------------------------#
# Linuxx8664gcc CMake System configuration file for the AliRoot Build System     #
#                                                                                # 
# Author: Johny Jose (johny.jose@cern.ch)                                        #
#                                                                                #
#--------------------------------------------------------------------------------#


cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#Retrieve Compiler Version

# Global Optimization
set(OPT)

# ------- Setting optimization flags for default configuration -------

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))
    set(DEFAULT_CXX_FLAGS "-O1 -g")
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

set(CXXWARN " ")

set(CXXFLAGS "${OPT} -fPIC -Dlinux")
set(CXXFLAGSNO "${NOOPT} -fPIC -Dlinux")
add_definitions(-Dlinux)

set(CFLAGS	"${OPT} -fPIC")
set(FFLAGS "${CLIBFOPT} ${FOPT} -fPIC")

set(SYSLIBS "-ldl -lcrypt -L/usr/X11R6/lib -lX11")

set(LDFLAGS "${OPT}")
set(SOFLAGS "${OPT} -shared")
set(ALLIB)



