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

# Global Optimization
set ( OPT )
# ------- Setting optimization flags for default configuration -------
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))
    set(DEFAULT_CXX_FLAGS "-O -g")
    set(OPT "${DEFAULT_CXX_FLAGS}")
    message("-- Setting compiler flags for default configuration: ${DEFAULT_CXX_FLAGS}")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "None"))
# --------------------------------------------------------------------
set ( NOOPT  "-g")

set ( CXXOPT  ${OPT})
set ( CXXNOOPT  ${NOOPT})
set ( COPT  ${OPT})
set ( FOPT  ${OPT})

set ( CLIBDEFS  "-DCERNLIB_SUN -DCERNLIB_BLDLIB -DCERNLIB_CZ")
set ( CLIBCXXOPTS )
set ( CLIBCOPT )
set ( CLIBFOPT  "${CLIBDEFS} -xpp=cpp")

set ( CXXFLAGS  "${CXXOPT} +w -KPIC -template=no%extdef")
set ( CXXFLAGSNO  "${CXXNOOPT} +w -KPIC -template=no%extdef")
set ( CFLAGS  "${COPT} -KPIC -erroff=%none")
set ( FFLAGS  "${FOPT} ${CLIBFOPT} ${CLIBDEFS}  -KPIC")
set ( DEPENDFFLAGS  ${FFLAGS})

set ( CINTFLAGS )

set ( LDFLAGS  "${OPT} -Qoption ld -t")

set ( SOFLAGS  "-G ${LDFLAGS} ${SHLIB}")

set ( SYSLIBS  "-L/usr/dt/lib -L/usr/openwin/lib -L/usr/ccs/lib -lXm -lXt -lX11 -lm -lgen -ldl -lsocket -lsunmath -lfsu -lfui -lnsl")

set ( EXEFLAGS  "-O -Qoption ld -t")
