#--------------------------------------------------------------------------------#
# Package File for                                                               #
# Author : Johny Jose (johny.jose@cern.ch)                                       #
# Variables Defined :                                                            #
#                                                                                #
# SRCS - C++ source files                                                        #
# HDRS - C++ header files                                                        #
# DHDR - ROOT Dictionary Linkdef header file                                     #
# CSRCS - C source files                                                         #
# CHDRS - C header files                                                         #
# EINCLUDE - Include directories                                                 #
# EDEFINE - Compiler definitions                                                 #
# ELIBS - Extra libraries to link                                                #
# ELIBSDIR - Extra library directories                                           #
# PACKFFLAGS - Fortran compiler flags for package                                #
# PACKCXXFLAGS - C++ compiler flags for package                                  #
# PACKCFLAGS - C compiler flags for package                                      #
# PACKSOFLAGS - Shared library linking flags                                     #
# PACKLDFLAGS - Module linker flags                                              #
# PACKBLIBS - Libraries to link (Executables only)                               #
# EXPORT - Header files to be exported                                           #
# CINTHDRS - Dictionary header files                                             #
# CINTAUTOLINK - Set automatic dictionary generation                             #
# ARLIBS - Archive Libraries and objects for linking (Executables only)          #
# SHLIBS - Shared Libraries and objects for linking (Executables only)           #
#--------------------------------------------------------------------------------#


execute_process (COMMAND ${CMAKE_C_COMPILER} -V 2>&1 | awk '{ if (NR == 1) print $$8}' | cut -d'.' -f1
                 OUTPUT_VARIABLE ICC_MAJOR
                 OUTPUT_STRIP_TRAILING_WHITESPACE) 
execute_process (COMMAND ${CMAKE_C_COMPILER} -V 2>&1 | awk '{ if(NR == 1) print $$8}' | cut -d'.' -f2
                 OUTPUT_VARIABLE ICC_MINOR
                 OUTPUT_STRIP_TRAILING_WHITESPACE)

# Global Optimization
set ( OPT  "-O -ip")
# ------- Setting optimization flags for default configuration -------
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))
    set(DEFAULT_CXX_FLAGS "-O -g")
    set(OPT "${DEFAULT_CXX_FLAGS}")
    message("-- Setting compiler flags for default configuration: ${DEFAULT_CXX_FLAGS}")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))
# --------------------------------------------------------------------
set ( NOOPT  "-O0")

set ( CXXOPT  "${OPT} -fPIC")
set ( CXXNOOPT  "${NOOPT} -fPIC")
set ( COPT  "${OPT} -fPIC")
set ( FOPT  "${OPT} -fPIC")

set ( CLIBDEFS  "-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ")
set ( CLIBCXXOPTS )
set ( CLIBCOPT )
set ( CLIBFOPT  ${CLIBDEFS})

set ( CXXFLAGS  ${CXXOPT})
set ( CXXFLAGSNO  ${CXXNOOPT})
set ( CFLAGS  ${COPT})
set ( FFLAGS  "${CLIBFOPT} ${FOPT}")
set ( DEPENDFFLAGS  ${FFLAGS})

set ( CINTFLAGS )

set ( LDFLAGS  ${OPT})

set ( SOFLAGS  "-shared")

set ( ALLIB )

execute_process(COMMAND find /lib -name 'libNoVersion*.so' | xargs --replace basename {} .so | sed -e 's/lib/ -l/'
                OUTPUT_VARIABLE LIBNOVER
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND which ifort | sed -e 's|bin|lib|' | sed -e 's|ifor||'
                OUTPUT_VARIABLE IFORT 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set (SYSLIBS "-lcrypt -L/usr/X11R6/lib -lX11 ${LIBNOVER} -L/usr/local/lib -lXt -L${IFORT} -lifcore -lifport")
                
