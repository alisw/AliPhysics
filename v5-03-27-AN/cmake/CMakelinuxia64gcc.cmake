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
set(NOOPT "-g -mieee -mno-soft-float")

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

set(CXXFLAGS "${OPT} -fPIC -pipe")
set(CXXFLAGSNO "${NOOPT} -fPIC -pipe")
set(CFLAGS	"${OPT} -Wall -fPIC -pipe -ansi")
set(FFLAGS "${CLIBFOPT} ${FOPT} -Wall -fPIC -pipe -fno-second-underscore")

set(SYSLIBS "-ldl -lg2c -lcrypt -L/usr/X11R6/lib -lX11")

set(LDFLAGS "${OPT}")
set(SOFLAGS "${OPT} -Wall -fPIC -pipe -shared -Wl")
set(SHLIB "-lg2c")
set(ALLIB)



