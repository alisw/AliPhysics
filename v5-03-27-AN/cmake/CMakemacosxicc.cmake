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

set(FINK_ROOT $ENV{FINK_ROOT})
if( NOT FINK_ROOT )
		
	set ( FINK_ROOT /usr/local)

endif( NOT FINK_ROOT )

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

set ( CLIBDEFS  "-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC")
set ( CLIBCXXOPTS )
set ( CLIBCOPT )
set ( CLIBFOPT  ${CLIBDEFS})

set ( CXXFLAGS  "${CXXOPT} -fPIC -wd1476 -I/usr/X11R6/include -I${FINK_ROOT}/include)"

set ( CXXFLAGSNO  ${CXXNOOPT})

set ( CFLAGS  "${COPT} -fPIC -restrict -I${FINK_ROOT}/include")

set ( FFLAGS  "${CLIBFOPT} ${FOPT}")

set ( CINTFLAGS )

set ( LDFLAGS  "${OPT} ${DICTLOAD}")

set ( SOFLAGS  "${OPT} -dynamiclib -undefined dynamic_lookup -single_module")

set ( DYFLAGS  "-dynamiclib -undefined dynamic_lookup -single_module")

set ( ALLIB )

set ( EXEFLAGS  "-bind_at_load")

execute_process(COMMAND find /lib -name 'libNoVersion*.so' | xargs --replace basename {} .so | sed -e 's/lib/ -l/'
                OUTPUT_VARIABLE LIBNOVER
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND which ifort | sed -e 's|bin|lib|'
                OUTPUT_VARIABLE IFORT 
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set (SYSLIBS "-lcrypt -L/usr/X11R6/lib -lX11 ${LIBNOVER} -L/usr/local/lib -lXt -L${IFORT} -lifcore -lifport")
