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

set ( CLIBDEFS  "-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ")
set ( CLIBCXXOPTS )
set ( CLIBCOPT )
set ( CLIBFOPT  ${CLIBDEFS})

set ( CXXFLAGS  "${OPT} -Wall -pipe -Woverloaded-virtual -Weffc++ -D_DLL")
set ( CXXFLAGSNO  "${NOOPT} -Wall -pipe -Woverloaded-virtual -Weffc++ -D_DLL")
set ( CFLAGS  "${OPT} -Wall -D_DLL")
set ( FFLAGS  "${CLIBFOPT} ${FOPT} -fno-second-underscore")

set ( CINTFLAGS )

set ( LDFLAGS  ${OPT})

set ( SOFLAGS  "${OPT} -shared -Wl,--export-all-symbols -Wl,-soname=$$TMPLIB -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc")

set ( SHLIB  ${ROOTLIBS} -lg2c)

set ( ALLIB )

set ( SYSLIBS  "-ldl -lg2c -lcrypt -L/usr/X11R6/lib -lX11")
