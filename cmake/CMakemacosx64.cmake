# -*- mode: cmake -*-

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
		
set ( OPT  "-O -g")
set ( NOOPT  "-O0 -g")

set ( CXXOPT  ${OPT})
set ( CXXNOOPT  ${NOOPT})
set ( COPT  ${OPT})

if( F77 MATCHES "g95")
	
	set ( FOPT "${FOPT}  -fbounds-check")

endif( F77 MATCHES "g95")

set ( CLIBDEFS  "-DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC")
set ( CLIBCXXOPTS )
set ( CLIBCOPT )
set ( CLIBFOPT  ${CLIBDEFS})

set ( CXXWARN  "-Wall -Wconversion -Wno-long-long -W -Weffc++ -Wshadow -Woverloaded-virtual -ansi")

set ( CXXSTF  "-pipe -fbounds-check -fsigned-char -fno-common -fmessage-length=0 -fno-default-inline -fno-inline -I/usr/X11R6/include -I${FINK_ROOT}/include")

set ( CXXFLAGS  "${CXXOPT} ${CXXSTF}")

set ( CXXFLAGSNO  "${CXXNOOPT} ${CXXSTF}")

set ( CFLAGS  "${COPT} -Wall -W -fno-common -pipe -I${FINK_ROOT}/include")

if(F77 MATCHES "g95")
	
	set ( FFLAGS "${FFLAGS}  -ftrace=full")
	set ( FFLAGS "${FFLAGS} -DFORTRAN_G95")

else()

	set ( FFLAGS "${FFLAGS} -DFORTRAN_GFORTRAN")

endif(F77 MATCHES "g95")

set ( CINTFLAGS )

set ( LDFLAGS  "${OPT} ${DICTLOAD}")

set ( SOFLAGS  "-m64 -dynamiclib -undefined dynamic_lookup -single_module")

set ( ALLIB )

set ( SYSLIBS  "-L/usr/X11R6/lib -lX11")


set ( EXEFLAGS  "-bind_at_load")

if( F77 MATCHES "g95")
	
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --print-search-dirs 
                  OUTPUT_VARIABLE SYSLIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX MATCH "^.*install:[^\n]*" SYSLIBS ${SYSLIBS})
  string(REGEX REPLACE "^.*install: " "" SYSLIBS ${SYSLIBS})
  set(SYSLIBS "-L${SHLIB} -lf95")
  set(DYLIB ${SYSLIBS})

else()

  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -m64 -print-file-name=libgfortran.dylib
                  OUTPUT_VARIABLE _shlib
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -m64 -print-file-name=libgfortranbegin.a
                  OUTPUT_VARIABLE _alib
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  set( SYSLIBS "${SYSLIBS} -ldl ${_shlib} ${_alib}")

endif(F77 MATCHES "g95")
