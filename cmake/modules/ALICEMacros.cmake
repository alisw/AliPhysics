# -*- mode: cmake -*-
#_______________________________________________________________________________
Macro(CHANGE_FILE_EXTENSION OLD_EXT NEW_EXT RESULT LIST)

# This is probably an obsolete Macro

  If (${OLD_EXT} MATCHES "^[*][.]+.*$")
    String(REGEX REPLACE "^[*]+([.].*)$" "\\1" OLD_EXT1 ${OLD_EXT}) 
  Endif (${OLD_EXT} MATCHES "^[*][.]+.*$")

  If (${NEW_EXT} MATCHES "^[*][.]+.*$")
    String(REGEX REPLACE "^[*]+([.].*)" "\\1" NEW_EXT1 ${NEW_EXT}) 
  Endif (${NEW_EXT} MATCHES "^[*][.]+.*$")

  Set(FLIST)
  Foreach (_current_FILE ${LIST})

    String(REGEX REPLACE "^(.*)${OLD_EXT1}$" "\\1${NEW_EXT1}" res ${_current_FILE})
    Set (FLIST ${FLIST} ${res})

  Endforeach (_current_FILE ${LIST})
  Set(${RESULT} ${FLIST})

Endmacro (CHANGE_FILE_EXTENSION)

#_______________________________________________________________________________
Macro (CHECK_OUT_OF_SOURCE_BUILD)

# Checks that the binary is built outside the source

   String(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
   If(insource)
      File(REMOVE_RECURSE ${CMAKE_SOURCE_DIR}/Testing)
      File(REMOVE ${CMAKE_SOURCE_DIR}/DartConfiguration.tcl)
      Message(FATAL_ERROR "ALIROOT should be installed as an out of source build, to keep the source directory clean. Please create a extra build directory and run the command 'cmake path_to_source_dir' in this newly created directory. You have also to delete the directory CMakeFiles and the file CMakeCache.txt in the source directory. Otherwise cmake will complain even if you run it from an out-of-source directory.") 
   Endif(insource)

EndMacro (CHECK_OUT_OF_SOURCE_BUILD)

#_______________________________________________________________________________
Function (AddLibrary LIB SRCS DHDRS)

# Adds an AliRoot library as a target

  Set(LDEF "${LIB}LinkDef.h")
  Set(DICT)
  If(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${LDEF})
# even with no cxx files, one may want to build an empty lib as a placeholder
# in AliRoot this is signalled by the existence of an (empty) ${LIB}LinkDef.h
    Set(DICT "G__${LIB}.cxx")
    Set(ASRCS ${SRCS} ${DICT})
    Root_Generate_Dictionary("${DHDRS}" "${LDEF}" "${DICT}" "${INCLUDE_DIRECTORIES}")
  Else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${LDEF})
    Message(STATUS "${LDEF} not found... probably building empty lib")
    Set(ASRCS ${SRCS})
  Endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${LDEF})

  Add_Library(${LIB} SHARED ${ASRCS})
  Target_Link_Libraries(${LIB} ${ALIROOT_LIBRARIES})
  Set_Target_Properties(${LIB} PROPERTIES ${ALIROOT_LIBRARY_PROPERTIES})
  
  Install(TARGETS ${LIB} DESTINATION ${ALIROOT_INSTALL_DIR}/lib
    COMPONENT shared)

  If(ALICE_STATIC_BUILD)
    Add_Library(${LIB}_a STATIC ${ASRCS})
    Install(TARGETS ${LIB}_a DESTINATION ${ALIROOT_INSTALL_DIR}/lib/static)
  EndIf(ALICE_STATIC_BUILD)

  If(ARGV3)
    Install(FILES ${ARGV3} DESTINATION ${ALIROOT_INSTALL_DIR}/include)
  Endif(ARGV3)

  CheckViols(${LIB} "${SRCS}")

EndFunction (AddLibrary)

#_______________________________________________________________________________
Macro (SetModule)

String(REGEX REPLACE "${ALICE_ROOT}/?([^/]*)/?$" "\\1" ALIROOT_MODULE "${CMAKE_CURRENT_SOURCE_DIR}")
Add_Definitions(-D_MODULE_=\"${ALIROOT_MODULE}\")

EndMacro(SetModule)


#_______________________________________________________________________________
Function (AddExecutable BIN SRCS LIBS)

# Adds an AliRoot executable as a target

  Add_Executable(${BIN} ${SRCS})
  Target_Link_Libraries(${BIN} ${ROOT_LIBRARIES} ${LIBS})
  Install(TARGETS ${BIN} DESTINATION ${ALIROOT_INSTALL_DIR}/bin)

  If(ALICE_STATIC_BUILD)
    Add_Executable(${BIN}_a ${SRCS})
    Set(_ar_libs)
    Foreach(_lib ${LIBS})
      Set(_ar_libs ${_ar_libs} ${_lib}_a)
    EndForeach(_lib ${LIBS})
    Foreach(_lib ${LIBS})
      Set(_ar_libs ${_ar_libs} ${_lib}_a)
    EndForeach(_lib ${LIBS})
    Target_Link_Libraries(${BIN}_a ${ROOT_LIBRARIES} ${_ar_libs})
    Install(TARGETS ${BIN}_a DESTINATION ${ALIROOT_INSTALL_DIR}/bin)
  EndIf(ALICE_STATIC_BUILD)

  CheckViols(${BIN} "${SRCS}")

EndFunction (AddExecutable)

#_______________________________________________________________________________
Macro (SetupSystem)

# Set up all system dependencies of AliRoot

Message(STATUS "Setting up system dependent parameters for ${ALICE_TARGET}" )

If(ALICE_TARGET STREQUAL macosx64) 

  Execute_process(
    COMMAND sw_vers -productVersion
    OUTPUT_VARIABLE MACOSX
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  String(REGEX REPLACE "^(.*)[.](.*)[.](.*)$" "\\1" MACOSX_MAJOR "${MACOSX}")
  String(REGEX REPLACE "^(.*)[.](.*)[.](.*)$" "\\2" MACOSX_MINOR "${MACOSX}")

  Message(STATUS "Mac OS X ${MACOSX_MAJOR}.${MACOSX_MINOR}")

  Find_Package(fink)
  
#  Set(CMAKE_CXX_COMPILER g++)
#  Set(CMAKE_C_COMPILER gcc)
# I found no way to make this work...
#  Set(CMAKE_CXX_LINK_EXECUTABLE 
#    "MACOSX_DEPLOYMENT_TARGET=${MACOSX_MAJOR}.${MACOSX_MINOR} ${CMAKE_CXX_LINK_EXECUTABLE}")

  Set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup -m64")

  Set(CMAKE_Fortran_FLAGS "-fno-second-underscore -m64")

  Set(CLIBDEFS "-DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC")

  Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -pipe -Wall -W -Wno-long-double -pipe -fbounds-check -fsigned-char -fno-common -fmessage-length=0 -Woverloaded-virtual -Weffc++ -Wconversion -Wshadow -fno-default-inline -fno-inline -I/usr/X11R6/include -I${FINK_ROOT}/include")
  Set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -m64 -Wall -W -fno-common -pipe -I${FINK_ROOT}/include")

  If(CMAKE_Fortran_COMPILER MATCHES g95) 
    Set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -ftrace=full -DFORTRAN_G95")
    Execute_process(COMMAND svn info $ENV{ALICE_ROOT} 
      COMMAND g95 --print-search-dirs 
      OUTPUT_VARIABLE _out)
    String(REGEX REPLACE "^.*install: *([^\n]*)/\n.*$" "\\1" _libdir ${_out})
    Set(ROOT_LIBRARIES "${ROOT_LIBRARIES} -L${_libdir} -lf95")
  Else(CMAKE_Fortran_COMPILER MATCHES g95)
    Set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DFORTRAN_GFORTRAN")
    Execute_process(
      COMMAND gfortran -m64 -print-file-name=libgfortran.dylib
      OUTPUT_VARIABLE FLIB
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    Set(ROOT_LIBRARIES "${ROOT_LIBRARIES} ${FLIB}")
    Execute_process(
      COMMAND gfortran -m64 -print-file-name=libgfortranbegin.a
      OUTPUT_VARIABLE FLIB
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    Set(ROOT_LIBRARIES "${ROOT_LIBRARIES} ${FLIB}")
  Endif(CMAKE_Fortran_COMPILER MATCHES g95) 

  Set(LINK_FLAGS "${LINK_FLAGS} -m64")

# I think this is useless
#  Set(ALIROOT_LIBRARIES "${ROOT_LIBRARIES} -L/usr/X11R6/lib -lX11")
  Set(ALIROOT_LIBRARIES "${ROOT_LIBRARIES}")

# Would like to use this, but did not manage on the Mac
#Include(FindOpenGL)
#Set(ROOT_LIBRARIES ${ROOT_LIBRARIES} ${OPENGL_LIBRARIES})
#Set(ROOT_INCLUDE_DIR ${ROOT_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR})
# Poor man's version of the above
  Set(ALIROOT_INCLUDE_DIR ${ROOT_INCLUDE_DIR} /usr/X11/include)

  Set(LINK_FLAGS "${LINK_FLAGS} -bind_at_load")

# 
# LD            = export MACOSX_DEPLOYMENT_TARGET=$(MACOSX_MAJOR).$(MACOSX_MINOR) ; \
# 		unset LD_PREBIND ; \
# 		g++
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 

Endif(ALICE_TARGET STREQUAL macosx64)

EndMacro (SetupSystem)

# ../build/Makefile.alphacxx6
# # -*- mode: makefile -*-
# # Makefile to build AliRoot for Alpha OSF1
# 
# # System dependent commands
# 
# XARGS = xargs
# 
# # The compilers
# CXX           = cxx 
# F77	      = f77
# 
# # Global optimisation
# OPT           = -O
# NOOPT	      = -O0 
# 
# CXXOPT	      = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_DECS -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = -I. $(CLIBDEFS)
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT)   -nostdnew -rtti -taso
# CXXFLAGSNO    = $(CXXNOOPT) -nostdnew -rtti -taso
# CFLAGS	      = $(COPT) -fPIC -pipe -mcpu=ev5 -D__osf__ -D__alpha
# FFLAGS        = $(FOPT) -nofor_main -warn noinformational -taso $(CLIBFOPT)
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(filter-out -warn noinformational,$(FFLAGS))
# 
# # rootcint flags
# CINTFLAGS     = -D__DECCXX
# 		
# LD            = cxx
# LDFLAGS       = 
# 
# SHLD	      = ld
# SOFLAGS       =  -L/usr/lib/cmplrs/cxx -rpath /usr/lib/cmplrs/cxx \
#                  -expect_unresolved "*" -msym -shared -taso \
#                  /usr/lib/cmplrs/cc/crt0.o /usr/lib/cmplrs/cxx/_main.o
# 
# SHLIB         = -lUfor -lfor -lFutil
# SOEXT 	      = so
# 
# #System libraries
# 
# # Flags for static libraries
# AFLAGS = $(filter-out -rpath /usr/lib/cmplrs/cxx -msym -shared /usr/lib/cmplrs/cc/crt0.o,$(SOFLAGS))
# AFLAGS += $(SHLIB)
# 
# # Additional flags and libraries for building aliroot executable
# SYSLIBS     := -lXm -lXt -lX11 -lPW -lUfor -lfor -lFutil -lots -taso -lbsd
# 
# # Cure funny problem 
# # sometimes in dependencies system include files of the sort
# # /usr/.../filename AND /usr/.../filename.cc are present
# # DEC believes that /usr/.../filename is the executable to be
# # built from /usr/.../filename.cc 
# # Just avoid this to happen
# 
# % : %.cc
# 	@;
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ../build/Makefile.hpuxacc
# # -*- mode: makefile -*-
# # Makefile to build AliRoot on HP-UX
# 
# # System dependent commands
# 
# XARGS = xargs
# 
# # The compilers
# CXX           = aCC
# F77	      = f77
# CC	      = cc
# 
# # Global optimisation
# OPT           = -g -O
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT          = $(OPT)
# FOPT          = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_HPUX -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      = +z -Ae 
# CLIBFOPT      = $(CLIBDEFS) -WF,-P
# 
# CXXFLAGS      = $(CXXOPT)   -Wc,-ansi_for_scope,on +Z -z +W70,495,740,749,823,829 -Dextname
# CXXFLAGSNO    = $(CXXNOOPT) -Wc,-ansi_for_scope,on +Z -z +W70,495,740,749,823,829 -Dextname
# CFLAGS	      = $(COPT) +Z -Ae
# FFLAGS        = $(CLIBFOPT) $(FOPT) +ppu +Z
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = aCC
# LDFLAGS       = $(OPT) -z -Wl,+s -Wl,-E
# 
# SHLD	      = $(LD)
# SOFLAGS       = -b /usr/lib/libm.sl 
# 
# SOEXT 	      = sl
# 
# # additional ROOT libraries
# SYSLIBS      := -lcrypt -L/usr/lib/X11R6 -lX11 
# 
# 
# 
# 
# 
# ../build/Makefile.linux
# # -*- mode: makefile -*-
# # Makefile to build AliRoot for Linux
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = g++ 
# CC	      = gcc
# CCMAJORV      = $(shell $(CC) -dumpversion | cut -d. -f1)
# CCMINORV      = $(shell $(CC) -dumpversion | cut -d. -f2)
# F77	      = $(shell root-config --f77)
# 
# # Global optimisation
# OPT           = -O -g
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# ifeq ($(CCMAJORV),2)
# CXXFLAGS       = $(OPT) -W -Wall -fPIC -pipe
# CXXFLAGSNO     = $(NOOPT) -W -Wall -fPIC -pipe
# else
# ifeq ($(CCMAJORV),3)
# CXXFLAGS       = $(OPT) -W -Wall -Weffc++ -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -Dlinux
# CXXFLAGSNO     = $(NOOPT) -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi
# else
# ifeq ($(CCMAJORV),4)
# CXXFLAGS       = $(OPT) -W -Wall -Weffc++ -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -Dlinux
# CXXFLAGSNO     = $(NOOPT) -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi
# else
# CXXFLAGS       = $(OPT) -W -Wall -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -ansi -Dlinux
# CXXFLAGSNO     = $(NOOPT) -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -ansi
# endif
# endif
# endif
# CFLAGS	       = $(OPT) -Wall -Werror -fPIC -pipe -Wno-long-long -pedantic-errors -ansi
# FFLAGS         = $(CLIBFOPT) $(FOPT) -fno-second-underscore
# 
# ifeq (g95,$(F77))
# FFLAGS	      +=-DFORTRAN_G95
# else
# ifeq (gfortran,$(F77))
# FFLAGS	      +=-DFORTRAN_GFORTRAN
# else
# FFLAGS	      +=
# endif
# endif
# 
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS   = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = g++
# LDFLAGS       = $(OPT) 
# 
# SHLD	      = $(LD)
# SOFLAGS       = $(OPT) -shared -Wl 
# SOEXT 	      = so
# 
# #System libraries
# SYSLIBS      := -ldl -lcrypt -L/usr/X11R6/lib -lX11
# 
# ifeq (g95,$(F77))
# SHLIB += -L$(shell g95 --print-search-dirs | sed -n -e 's/install: //p') -lf95
# else
# ifeq (gfortran,$(F77))
# SHLIB := $(shell gfortran -print-file-name=libgfortran.so)
# SHLIB += $(shell gfortran -print-file-name=libgfortranbegin.a)
# SYSLIBS += $(SHLIB)
# else
# SHLIB         = -lg2c
# SYSLIBS +=  -lg2c
# endif
# endif
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# # additional ROOT libraries
# 
# 
# 
# 
# 
# 
# ../build/Makefile.linuxalphagcc
# # -*- mode: makefile -*-
# # Makefile to build AliRoot for Linux on alpha
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = g++ 
# F77	      = g77
# CC	      = gcc
# CCMAJORV      = $(shell $(CC) -dumpversion | cut -d. -f1)
# CCMINORV      = $(shell $(CC) -dumpversion | cut -d. -f2)
# 
# # Global optimisation
# OPT           = -O -g -mieee -mno-soft-float
# NOOPT         = -g -mieee -mno-soft-float
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_UNIX -DCERNLIB_DECS -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# 
# CXXFLAGS       = $(OPT) -Wall -fPIC -pipe
# CXXFLAGSNO     = $(NOOPT) -Wall -fPIC -pipe
# CFLAGS	       = $(OPT) -Wall -fPIC -pipe -ansi
# FFLAGS         = $(CLIBFOPT) $(FOPT) -Wall -fPIC -pipe -fno-second-underscore
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS   = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = g++
# LDFLAGS       = $(OPT) 
# 
# SHLD	      = $(LD)
# SOFLAGS       = $(OPT) -Wall -fPIC -pipe -shared -Wl
# SHLIB         = -lg2c
# SOEXT 	      = so
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# #System libraries
# SYSLIBS      := -ldl -lg2c -lcrypt -L/usr/X11R6/lib -lX11
# 
# 
# 
# 
# 
# ../build/Makefile.linuxia64ecc
# # -*- mode: makefile -*-
# # Makefile for AliRoot for Itanium/Linux with Intel icc compiler
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = icc 
# F77	      = ifort
# CC	      = icc
# 
# # Global optimisation
# OPT           = -g #-O
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(OPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT)
# CXXFLAGSNO    = $(CXXNOOPT)
# CFLAGS	      = $(COPT)
# FFLAGS        = $(CLIBFOPT) $(FOPT)
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = icpc
# LDFLAGS       = $(OPT) 
# 
# SHLD	      = $(LD)
# SOFLAGS       = -Wl,-soname,$$TMPLIB -shared -O
# SHLIB         =
# SOEXT 	      = so
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# # additional ROOT libraries
# 
# LIBNOVER      = `find /lib -name 'libNoVersion*.so' | xargs --replace basename {} .so | sed -e 's/lib/ -l/'`
# 
# SYSLIBS      := $(LIBNOVER) -lcrypt -L/usr/local/lib -lXt -lCEPCF90 -lF90 \
#                 -lPEPCF90 -lintrins -L/usr/X11R6/lib -lX11
# ../build/Makefile.linuxia64gcc
# # -*- mode: makefile -*-
# # Makefile for AliRoot for Itanium/Linux with gcc
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = g++ 
# F77	      = g77
# CC	      = gcc
# 
# # Global optimisation
# OPT           = -g -O
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT)  -Wall -Weffc++ -fPIC -pipe
# CXXFLAGSNO    = $(CXXNOOPT)  -Wall -Weffc++ -fPIC -pipe
# CFLAGS	      = -Wall -fPIC -pipe
# FFLAGS        = $(CLIBFOPT) $(FOPT) -fno-second-underscore
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = g++
# LDFLAGS       = $(OPT) -Wl,-Map -Wl,$@.map
# 
# SHLD	      = $(LD)
# SOFLAGS       = -shared -Wl -O2
# SOEXT 	      = so
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# # additional ROOT libraries
# 
# LIBNOVER      = `find /lib -name 'libNoVersion*.so' | xargs --replace basename {} .so | sed -e 's/lib/ -l/'`
# 
# #System libraries
# SYSLIBS      := -ldl -lbsd -lg2c -L/usr/X11R6/lib -lX11  $(LIBNOVER)
# 
# 
# 
# 
# 
# ../build/Makefile.linuxicc
# # -*- mode: makefile -*-
# # Makefile for AliRoot for Linux with the Intel icc compiler
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = icc
# CC            = icc
# 
# # Compiler version:
# ICC_MAJOR    := $(shell $(CXX) -V 2>&1 | awk '{ if (NR==1) print $$8 }' | \
#                 cut -d'.' -f1)
# ICC_MINOR    := $(shell $(CXX) -V 2>&1 | awk '{ if (NR==1) print $$8 }' | \
#                 cut -d'.' -f2)
# 
# F77           = ifort
# 
# # Global optimisation
# OPT           = -O -ip
# NOOPT         = -O0
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT          = $(OPT)
# FOPT          = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT)
# CXXFLAGSNO    = $(CXXNOOPT)
# CFLAGS        = $(COPT)
# FFLAGS        = $(CLIBFOPT) $(FOPT)
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     =
# 
# LD            = icpc
# LDFLAGS       = $(OPT)
# 
# SHLD          = $(LD)
# SOFLAGS       = -Wl,-soname,$$TMPLIB -shared $(OPT)
# SHLIB         =
# SOEXT 	      = so
# 
# ALLD          = ar
# ALFLAGS       = cr
# ALLIB         =
# AEXT 	      = a
# 
# # additional ROOT libraries
# 
# LIBNOVER      = `find /lib -name 'libNoVersion*.so' | xargs --replace basename {} .so | sed -e 's/lib/ -l/'`
# 
# #System libraries
# SYSLIBS      := -lcrypt -L/usr/X11R6/lib -lX11 $(LIBNOVER) -L/usr/local/lib \
#                 -lXt -L$(shell which ifort | sed -e 's|bin/ifort|lib|') \
#                 -lifcore -lifport
# ../build/Makefile.linuxx8664gcc
# # -*- mode: makefile -*-
# # Makefile to build AliRoot for Linux
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = g++ 
# CC	      = gcc
# CCMAJORV      = $(shell $(CC) -dumpversion | cut -d. -f1)
# CCMINORV      = $(shell $(CC) -dumpversion | cut -d. -f2)
# F77	      = $(shell root-config --f77)
# 
# # Global optimisation
# OPT           = -O0 -g
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# ifeq ($(CCMAJORV),2)
# CXXFLAGS       = $(OPT) -W -Wall -fPIC -pipe -m64
# CXXFLAGSNO     = $(NOOPT) -W -Wall -fPIC -pipe -m64
# else
# ifeq ($(CCMAJORV),3)
# CXXFLAGS       = $(OPT) -W -Wall -Weffc++ -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -Dlinux -m64
# CXXFLAGSNO     = $(NOOPT) -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -m64
# else
# ifeq ($(CCMAJORV),4)
# CXXFLAGS       = $(OPT) -W -Wall -Weffc++ -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -Dlinux -m64
# CXXFLAGSNO     = $(NOOPT) -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -m64
# else
# CXXFLAGS       = $(OPT) -W -Wall -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -ansi -Dlinux -m64
# CXXFLAGSNO     = $(NOOPT) -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -ansi -m64
# endif
# endif
# endif
# CFLAGS	       = $(OPT) -Wall -Werror -fPIC -pipe -Wno-long-long -pedantic-errors -ansi -m64 
# FFLAGS         = $(CLIBFOPT) $(FOPT) -fno-second-underscore -fPIC -fno-f2c -m64
# 
# ifeq (g95,$(F77))
# FFLAGS	      +=-DFORTRAN_G95
# else
# ifeq (gfortran,$(F77))
# FFLAGS	      +=-DFORTRAN_GFORTRAN
# else
# FFLAGS	      +=
# endif
# endif
# 
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS   = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = g++
# LDFLAGS       = $(OPT) 
# 
# SHLD	      = $(LD) 
# SOFLAGS       = $(OPT) -shared -Wl 
# SOEXT 	      = so
# 
# #System libraries
# LIBNOVER      = `find /lib64 -name 'libNoVersion*.so' | xargs --replace basename {} .so | sed -e 's/lib64/ -l/'`
# 
# SYSLIBS      := -ldl -lcrypt -L/usr/X11R6/lib64 -lX11  $(LIBNOVER)
# 
# ifeq (g95,$(F77))
# SHLIB += -L$(shell g95 --print-search-dirs | sed -n -e 's/install: //p') -lf95
# else
# ifeq (gfortran,$(F77))
# SHLIB := $(shell gfortran -print-file-name=libgfortran.so)
# SHLIB += $(shell gfortran -print-file-name=libgfortranbegin.a)
# else
# SHLIB         = -lg2c
# SYSLIBS +=  -lg2c
# endif
# endif
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ../build/Makefile.macosx
# # -*- mode: makefile -*-
# # Makefile for AliRoot for MacOS X with gcc
# 
# XARGS = xargs
# 
# # OS version
# MACOSX_MAJOR := $(strip $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 1))
# MACOSX_MINOR := $(strip $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2))
# 
# # fink directories
# FINK_ROOT := $(shell which fink | sed -e 's?/bin/fink??')
# ifeq (,$(FINK_ROOT))
# # No fink, build will probably fail, but we try a guess
# FINK_ROOT=/usr/local
# endif
# 
# # The compilers
# CXX           = g++ 
# CC	      = gcc
# F77	      = $(shell root-config --f77)
# 
# # Global optimisation
# OPT           = -O -g
# NOOPT         = -O0 -g
# 
# CXXOPT        = $(OPT) 
# CXXNOOPT      = $(NOOPT) 
# COPT	      = $(OPT)
# FOPT	      = $(OPT) -fno-second-underscore 
# ifeq (g95,$(F77))
# FOPT	     += -fbounds-check
# endif
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# CXXSTF        =  -pipe -Wall -W -Wno-long-double -pipe -fbounds-check -fsigned-char -fno-common -fmessage-length=0 -Woverloaded-virtual -Weffc++ -Wconversion -Wshadow -fno-default-inline -fno-inline -I/usr/X11R6/include -I$(FINK_ROOT)/include
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT) $(CXXSTF)
# 
# CXXFLAGSNO    = $(CXXNOOPT) $(CXXSTF) 
# 
# CFLAGS	      = $(COPT) -Wall -W -fno-common -pipe -I$(FINK_ROOT)/include
# 
# FFLAGS        = $(CLIBFOPT) $(FOPT)
# ifeq (g95,$(F77))
# FFLAGS        += -ftrace=full
# FFLAGS        +=-DFORTRAN_G95
# else
# FFLAGS        +=-DFORTRAN_GFORTRAN
# endif
# 
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = export MACOSX_DEPLOYMENT_TARGET=$(MACOSX_MAJOR).$(MACOSX_MINOR) ; \
# 		unset LD_PREBIND ; \
# 		g++
# LDFLAGS       = $(OPT) $(DICTLOAD)
# 
# SHLD	     := $(LD)
# SOFLAGS      := -bundle -undefined dynamic_lookup
# SHLIB        := 
# SOEXT 	     := so
# 
# DYLD	     := $(LD)
# DYFLAGS       = -dynamiclib -undefined dynamic_lookup -single_module
# DYLIB        := 
# DYEXT        := dylib
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# DEPENDCXXFLAGS = $(CXXFLAGS) -I/usr/include/sys
# 
# SYSLIBS      := -L/usr/X11R6/lib -lX11
# 
# EXEFLAGS     := -bind_at_load 
# 
# #System libraries
# 
# ifeq (g95,$(F77))
# SYSLIBS += -L$(shell g95 --print-search-dirs | sed -n -e 's/install: //p') -lf95
# DYLIB += -L$(shell g95 --print-search-dirs | sed -n -e 's/install: //p') -lf95
# else
# SYSLIBS += -ldl $(shell $(F77) -print-file-name=libgfortran.dylib)
# SYSLIBS += $(shell $(F77) -print-file-name=libgfortranbegin.a)
# endif
# ../build/Makefile.macosxicc
# # -*- mode: makefile -*-
# # Makefile for AliRoot for MacOS X with gcc
# 
# XARGS = xargs
# 
# # fink directories
# FINK_ROOT := $(shell which fink | sed -e 's?/bin/fink??')
# ifeq (,$(FINK_ROOT))
# # No fink, build will probably fail, but we try a guess
# FINK_ROOT=/usr/local
# endif
# 
# # The compilers
# CXX           = icc 
# CC	      = icc
# 
# F77	      = ifort
# 
# # Global optimisation
# OPT           = -O -g
# NOOPT         = -g
# 
# CXXOPT        = $(OPT) 
# CXXNOOPT      = $(NOOPT) 
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT) -fPIC -wd1476 -I/usr/X11R6/include -I$(FINK_ROOT)/include 
# 
# CXXFLAGSNO    = $(CXXNOOPT) 
# 
# CFLAGS	      = $(COPT) -fPIC -restrict -I$(FINK_ROOT)/include
# 
# FFLAGS        = $(CLIBFOPT) $(FOPT)
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = icpc
# LDFLAGS       = $(OPT) $(DICTLOAD)
# 
# SHLD	     := $(LD)
# SOFLAGS      := $(OPT) -dynamiclib -undefined dynamic_lookup -single_module
# SHLIB        := 
# SOEXT 	     := so
# 
# DYLD	     := $(LD)
# DYFLAGS       = -dynamiclib -undefined dynamic_lookup -single_module
# DYLIB        := 
# DYEXT        := dylib
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# DEPENDCXXFLAGS = $(CXXFLAGS) -I/usr/include/sys
# 
# EXEFLAGS     := -bind_at_load 
# 
# #System libraries
# 
# SYSLIBS      := -L/usr/X11R6/lib -lX11 -lGL $(LIBNOVER) -L/usr/local/lib \
#                 -lXt -L$(shell which ifort | sed -e 's|bin/ifort|lib|') \
#                 -lifcore -lifport
# 
# 
# ../build/Makefile.macosxxlc
# # -*- mode: makefile -*-
# # Makefile for AliRoot for MacOS X using the IBM xlc compiler
# 
# # The compilers
# CXX           = xlC
# F77	      = xlf
# CC	      = xlc
# 
# # Global optimisation
# OPT           = -qnoopt #-O3 -g
# NOOPT         = -qnoopt
# 
# CXXOPT        = $(OPT) 
# CXXNOOPT      = $(NOOPT) 
# COPT	      = $(OPT)
# FOPT	      = $(OPT) 
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(shell echo $(CLIBDEFS) | sed -e's/-D/-WF,-D/g')
# 
# # Compiler flags
# CXXFLAGS      = $(CXXOPT) -qpic -I/sw/include \
#                 -qflttrap=overflow:zerodivide:invalid:inexact:enable
# CXXFLAGSNO    = $(CXXNOOPT) -qpic -I/sw/include \
#                 -qflttrap=overflow:zerodivide:invalid:inexact:enable
# CFLAGS	      = -qpic -qextname -I/sw/include \
#                 -qflttrap=overflow:zerodivide:invalid:inexact:enable
# FFLAGS        = $(CLIBFOPT)  $(FOPT) -qpic \
#                 -qflttrap=overflow:zerodivide:invalid:inexact:enable
# # rmkdepend flags for building dependencies 
# DEPENDFFLAGS  = $(FFLAGS)
# DEPENDCXXFLAGS = $(CXXFLAGS) -I/usr/include/sys
# 
# # rootcint flags
# CINTFLAGS     = 
# 		
# LD            = xlC
# LDFLAGS       = $(OPT) 
# 
# SHLD	     := export MACOSX_DEPLOYMENT_TARGET=10.3 ; \
# 		unset LD_PREBIND ; \
# 		xlC
# SOFLAGS      := -bundle -undefined dynamic_lookup
# #SHLIB        := -lg2c
# SHLIB        :=
# SOEXT 	     := so
# 
# DYLD	     := export MACOSX_DEPLOYMENT_TARGET=10.3 ; \
# 		unset LD_PREBIND ; \
# 		xlC
# DYFLAGS       = -qmkshrobj -undefined dynamic_lookup -single_module
# DYLIB        :=
# DYEXT        := dylib
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# 
# #System libraries
# SYSLIBS      := -ldl -L/usr/X11R6/lib -lX11
# 
# EXEFLAGS     := -Wl,-bind_at_load
# 
# 
# 
# 
# ../build/Makefile.solarisCC5
# # -*- mode: makefile -*-
# # Makefile to build ALIROOT for SunOS
# 
# # System dependent commands
# 
# XARGS = xargs
# 
# # The compilers
# CXX           = /opt/SUNWspro/bin/CC
# CC	      = /opt/SUNWspro/bin/cc
# F77	      = /opt/SUNWspro/bin/f77
# 
# # Global optimisation
# OPT           = -g -O 
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_SUN -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS) -xpp=cpp
# 		
# # Compiler flags
# CXXFLAGS      = $(CXXOPT) +w -KPIC -features=rtti -D_XOPEN_SOURCE -D_XOPEN_VERSION=4 -D__EXTENSIONS__
# CXXFLAGSNO    = $(CXXNOOPT) +w -KPIC -features=rtti
# CFLAGS	      = $(COPT) -KPIC -erroff=%none
# FFLAGS        = $(FOPT) $(CLIBFOPT) $(CLIBDEFS)  -KPIC
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS  = $(FFLAGS) 
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# SHLIB         = 
# SOEXT 	      = so
# 
# LD            = /opt/SUNWspro/bin/CC
# LDFLAGS       = $(OPT) -Qoption ld -t
# 
# SHLD	      = $(LD)
# SOFLAGS       = -G $(LDFLAGS) $(SHLIB)  
# 
# SYSLIBS      := -L/usr/dt/lib -L/usr/openwin/lib -L/usr/ccs/lib -lXm -lXt \
#                 -lX11 -lm -lgen -ldl -lsocket -lsunmath -lfsu -lfui -lnsl
# 
# # Additional flags and libraries for building aliroot executable
# 
# EXEFLAGS     := -O -Qoption ld -t 
# 
# 
# 
# 
# 
# ../build/Makefile.win32gcc
# # -*- mode: makefile -*-
# # Makefile to build AliRoot for Linux
# 
# # System dependent commands
# 
# XARGS = xargs -r
# 
# # The compilers
# CXX           = g++ 
# F77	      = g77
# CC	      = gcc
# CCMAJORV      = $(shell $(CC) -dumpversion | cut -d. -f1)
# CCMINORV      = $(shell $(CC) -dumpversion | cut -d. -f2)
# 
# # Global optimisation
# OPT           = -O -g
# NOOPT         = -g
# 
# CXXOPT        = $(OPT)
# CXXNOOPT      = $(NOOPT)
# COPT	      = $(OPT)
# FOPT	      = $(OPT)
# 
# # CERNLIB defines
# CLIBDEFS      = -DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ
# CLIBCXXOPTS   =
# CLIBCOPT      =
# CLIBFOPT      = $(CLIBDEFS)
# 
# # Compiler flags
# CXXFLAGS       = $(OPT) -Wall -pipe -Woverloaded-virtual -Weffc++ -D_DLL
# CXXFLAGSNO     = $(NOOPT) -Wall -pipe -Woverloaded-virtual -Weffc++ -D_DLL
# CFLAGS	       = $(OPT) -Wall -D_DLL
# FFLAGS         = $(CLIBFOPT) $(FOPT) -fno-second-underscore
# # rmkdepend flags for building dependencies of FORTRAN files
# DEPENDFFLAGS   = $(FFLAGS)
# 
# # rootcint flags
# CINTFLAGS     = 
# 
# LD            = g++
# LDFLAGS       = $(OPT) 
# 
# SHLD	      = $(LD)
# SOFLAGS       = $(OPT) -shared -Wl,--export-all-symbols -Wl,-soname=$$TMPLIB -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc
# SHLIB         = $(shell root-config --libs) -lg2c
# SOEXT 	      = dll
# 
# ALLD	      = ar
# ALFLAGS       = cr
# ALLIB         = 
# AEXT 	      = a
# 
# # additional ROOT libraries
# 
# #System libraries
# SYSLIBS      := -ldl -lg2c -lcrypt -L/usr/X11R6/lib -lX11
# 
# 
# 
# 
# 






