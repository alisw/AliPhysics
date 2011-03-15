#--------------------------------------------------------------------------------#
# Maxosx CMake System configuration file for the AliRoot Build System            #
#                                                                                #
# Author : Johny Jose (johny.jose@cern.ch)                                       #
#                                                                                #
#--------------------------------------------------------------------------------#

execute_process (COMMAND sw_vers | sed -n 's/ProductVersion://p' | cut -d. -f1
                 OUTPUT_VARIABLE MACOSX_MAJOR
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process (COMMAND sw_vers | sed -n 's/ProductVersion://p' | cut -d. -f2
                 OUTPUT_VARIABLE MACOSX_MINOR
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

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

set ( CXXOPT  "${OPT}")
set ( CXXNOOPT  "${NOOPT}")
set ( COPT  "${OPT}")

set ( CLIBDEFS  "-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC")
set ( CLIBCXXOPTS )
set ( CLIBCOPT )
set ( CLIBFOPT  ${CLIBDEFS})

set ( CXXWARN  "-Wall -Wno-long-long -W -Weffc++ -Wshadow -Woverloaded-virtual -ansi")

set ( CXXSTF  "-pipe -fbounds-check -fsigned-char -fno-common -fmessage-length=0 -fno-default-inline -fno-inline -I/usr/X11R6/include -I${FINK_ROOT}/include")

set ( CXXFLAGS  "${CXXOPT} ${CXXSTF}")

set ( CXXFLAGSNO  "${CXXNOOPT} ${CXXSTF}")

set ( CFLAGS  "${COPT} -Wall -W -fno-common -pipe -I${FINK_ROOT}/include")

set ( FFLAGS "${FFLAGS} -DFORTRAN_GFORTRAN")

set ( DEPENDFFLAGS  ${FFLAGS})

set ( CINTFLAGS )

set ( LDFLAGS  "${OPT} ${DICTLOAD}")

set ( SOFLAGS  "-dynamiclib -undefined dynamic_lookup -single_module")

set ( ALLIB )

set ( DEPENDCXXFLAGS  "${CXXFLAGS} -I/usr/include/sys")

set ( SYSLIBS  "-L/usr/X11R6/lib -lX11")

set ( EXEFLAGS  "-bind_at_load")

set (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")

set (CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem")

 
execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.dylib
                  OUTPUT_VARIABLE _shlib
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortranbegin.a
                  OUTPUT_VARIABLE _alib
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
set( SYSLIBS "${SYSLIBS} -ldl ${_shlib} ${_alib}")

