dnl -*- mode: autoconf -*- 
dnl
dnl $Id$
dnl template for the configuration script for the Alice HLT 
dnl analysis tutorial and example library
dnl 
dnl ------------------------------------------------------------------
AC_INIT([Alice High Level Trigger Tutorial] , [0.1], [Matthias.Richter@ift.uib.no], alice-hlt-sample)

dnl ------------------------------------------------------------------
AC_CANONICAL_SYSTEM
AC_PREFIX_DEFAULT(${PWD})
AC_CONFIG_SRCDIR(AliHLTDummyComponent.cxx)
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_DEBUG
AC_OPTIMIZATION
AM_CONDITIONAL(STANDALONE_SAMPLELIB, test 1)

dnl ------------------------------------------------------------------
dnl
dnl Check for ROOT
dnl
ROOT_PATH(, [have_root=1], [AC_ERROR([Stop! The HLT package needs ROOT.])])
AM_CONDITIONAL(HAVE_ROOT, test "x$have_root" = "x1" )
AC_SUBST([ROOTSYS])

dnl test for additional required root libraries
  # we check for the libSTEER library which is linked to
  # - the ROOT libs libGeom libMinuit libVMC libEG
  # - the AliRoot libESD libRAWData (libRAWDatarec from v4-04-Rev-07)
  # - ROOT libCint needs also libdl
  # - from Nov 1 2006 TTreeFormula is needed by AliTagAnalysis and requires
  #   libTreePlayer.so
  # - from Jan 07 libESD also depends on libXMLIO
LIBS='-ldl'
if test "x$have_root" = "x1"; then
  save_CPPFLAGS=$CPPFLAGS
  save_LDFLAGS=$LDFLAGS
  save_LIBS=$LIBS
  
  ROOT_CHECKLIBS='Geom Minuit EG VMC TreePlayer XMLIO'
  for CHECKLIB in $ROOT_CHECKLIBS ; do
    CPPFLAGS="$save_CPPFLAGS"
    LDFLAGS="$save_LDFLAGS -L${ROOTLIBDIR}"
    LIBS="$save_LIBS $ROOTLIBS $ADD_ROOTLIBS"
    AC_CHECK_LIB([$CHECKLIB],[_init], [ADD_ROOTLIBS="$ADD_ROOTLIBS -l$CHECKLIB"])
  done
  
  CPPFLAGS=$save_CPPFLAGS
  LDFLAGS=$save_LDFLAGS
  LIBS=$save_LIBS
fi

CPPFLAGS="$CPPFLAGS ${ROOTCFLAGS}"

dnl ------------------------------------------------------------------
dnl check for AliRoot features
AC_LANG_PUSH(C++)
have_aliroot=no
AC_ARG_WITH(aliroot,[  --with-aliroot	  top of the AliRoot installation],
		    [test -d $with_aliroot && ALICE_ROOT=$with_aliroot],
		    [])

if test "x$ALICE_ROOT" != "x" && test -d $ALICE_ROOT ; then
  dnl ------------------------------------------------------------------
  dnl Try to estimate the system architecture
  case $host_os:$host_cpu in
  linux*:x86_64*)	alice_target='linuxx8664gcc'		;;
  linux*)		alice_target='linux'			;;
  *)			alice_target='unknown'			;;
  esac
  if test "x$alice_target" = "xunknown" ; then
    if test -z $ALICE_TARGET ; then
    AC_MSG_ERROR([Can not estimate system architecture.
    To avoid the problem, set the ALICE_TARGET variable appropriately.
    Please send the following information to Matthias.Richter@ift.uib.no:
	${PACKAGE}: no AliRoot target definition for os=$host_os and cpu=$host_cpu ])
    else
    AC_MSG_NOTICE([Unknown system architecture.
    Please send the following information to Matthias.Richter@ift.uib.no:
	${PACKAGE}: no AliRoot target definition for os=$host_os and cpu=$host_cpu ])  
    fi
  else
    if test "x$ALICE_TARGET" != "x" && test "x$ALICE_TARGET" != "x$alice_target" ; then
    AC_MSG_WARN([The ALICE_TARGET variable does not match your current system
    overriding $ALICE_TARGET by $alice_target])
    fi
    ALICE_TARGET=$alice_target
  fi
  ALIROOTLIBDIR=${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}
  ALIROOTINCDIR=${ALICE_ROOT}/include
  test -d ${ALIROOTLIBDIR} || AC_MSG_WARN([can not find AliRoot library directory $ALIROOTLIBDIR])
  test -d ${ALIROOTINCDIR} || AC_MSG_WARN([can not find AliRoot include directory $ALIROOTINCDIR])
fi

AC_MSG_CHECKING([for AliRoot])
if test "x$ALICE_ROOT" != "x" \
   && test -d ${ALIROOTLIBDIR} \
   && test -d ${ALIROOTINCDIR}; then
  have_aliroot=$ALICE_ROOT
else
  ALIROOTLIBDIR=
  ALIROOTINCDIR=
fi
AC_MSG_RESULT([$have_aliroot])

if test ! "x$have_aliroot" = "xno" ; then
  ALIROOT_CPPFLAGS="-I${ALIROOTINCDIR} -I${ALICE_ROOT}/RAW"
  ALIROOT_LDFLAGS="-L${ALIROOTLIBDIR}"
  ALIROOT_LIBS="-lESD $ADD_ROOTLIBS"
  save_CPPFLAGS=$CPPFLAGS
  save_LDFLAGS=$LDFLAGS
  save_LIBS=$LIBS
  CPPFLAGS="$save_CPPFLAGS $ALIROOT_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS -L${ROOTLIBDIR} $ALIROOT_LDFLAGS"
  
  # check for certain AliRoot libraries/files/features
  # CBD library is present since AliRoot version v4-05-00 (02.06.2006)
  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  CHECKLIB=CDB
  AC_CHECK_LIB([$CHECKLIB],[_init],[ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"])

  # splitted RAW libraries since AliRoot version v4-04-Rev-07 (09.08.2006)
  have_alirawdata=no
  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  CHECKLIB=RAWDatabase
  AC_CHECK_LIB([$CHECKLIB],[_init],
	[ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"
	 LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
	 CHECKLIB=RAWDatarec
	 AC_CHECK_LIB([$CHECKLIB],[_init],
               [ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"
	        AC_DEFINE(HAVE_ALIRAWDATA)
	        have_alirawdata=$CHECKLIB])],
	[# second pass with -RAWData
	 CHECKLIB="-lRAWData"
	 LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS $CHECKLIB -lSTEER"
	 AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <AliRawReaderMemory.h>
 				       #include <AliRawReaderFile.h>
 				       #include <AliRawReaderDate.h>
 				       #include <AliRawReaderRoot.h>], 
                                     [AliRawReaderMemory mreader;
 				       AliRawReaderFile freader;
 				       AliRawReaderDate dreader(NULL,0);
 				       AliRawReaderRoot rreader(NULL,0);])],
                                     [AC_DEFINE(HAVE_ALIRAWDATA)
 				       ALIROOT_LIBS="$ALIROOT_LIBS $CHECKLIB"
 				       have_alirawdata=$CHECKLIB], 
 				      [have_alirawdata=no])
	]) dnl AC_CHECK_LIB RAWDatabase
  AC_MSG_CHECKING(for AliRawReader classes in RAWData libraries)
  AC_MSG_RESULT($have_alirawdata)

  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  CHECKLIB=STEER
  AC_CHECK_LIB([$CHECKLIB],[_init],[ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"])

  if test "x$have_aliroot" = "xno" ; then
    ALIROOT_CPPFLAGS=
    ALIROOT_LDFLAGS=
    ALIROOT_LIBS=
  fi
  CPPFLAGS=$save_CPPFLAGS
  LDFLAGS=$save_LDFLAGS
  LIBS=$save_LIBS
fi

dnl ------------------------------------------------------------------
dnl check for the HLT component framework
AC_MSG_CHECKING([for HLT component framework in ])
AC_ARG_WITH(hltbase, [  --with-hltbase	  installation path of the HLT component framework],
		     [],
		     [test -n $have_aliroot && with_hltbase=$have_aliroot ])
if test -z $with_hltbase || ! test -d $with_hltbase ; then
  with_hltbase=no
fi
AC_MSG_RESULT([$with_hltbase])
if test "x$with_hltbase" = "xno" || test "x$with_hltbase" = "x$ALICE_ROOT" ; then
  HLTBASE_CPPFLAGS=-I${ALICE_ROOT}/HLT/BASE
  HLTBASE_LDFLAGS=
else
  HLTBASE_CPPFLAGS=-I${with_hltbase}/include/alice-hlt
  HLTBASE_LDFLAGS=-L${with_hltbase}/lib
fi

  save_CPPFLAGS=$CPPFLAGS
  save_LDFLAGS=$LDFLAGS
  save_LIBS=$LIBS
  CPPFLAGS="$save_CPPFLAGS $HLTBASE_CPPFLAGS -L${ALIROOTLIBDIR}"
  LDFLAGS="$save_LDFLAGS -L${ROOTLIBDIR} $ALIROOT_LDFLAGS $HLTBASE_LDFLAGS"
  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  AC_CHECK_HEADER([AliHLTComponent.h])
  AC_CHECK_HEADER([AliHLTLogging.h])
  CHECKLIB=HLTbase
  AC_CHECK_LIB([$CHECKLIB],[_init],
	       [have_hltbase=yes
	        HLTBASE_LDFLAGS="$HLTBASE_LDFLAGS -l$CHECKLIB"])

  test "x$have_hltbase" = "xyes" || \
     AC_ERROR([Stop! The package needs the HLT component framework.
      Use  --with-hltbase option to point to the HLT component framework
      installation])

  CPPFLAGS=$save_CPPFLAGS
  LDFLAGS=$save_LDFLAGS
  LIBS=$save_LIBS

AC_LANG_POP(C++)
AC_SUBST([ALIROOT_CPPFLAGS])
AC_SUBST([ALIROOT_LDFLAGS])
AC_SUBST([ALIROOTLIBDIR])
AC_SUBST([ALIROOT_LIBS])
AC_SUBST([HLTBASE_CPPFLAGS])
AC_SUBST([HLTBASE_LDFLAGS])

dnl ------------------------------------------------------------------
dnl
dnl Documentation
dnl
AC_ARG_VAR(DOXYGEN, The Documentation Generator)
AC_PATH_PROG(PERL, perl)
AC_PATH_PROG(DOXYGEN, doxygen)
AM_CONDITIONAL(HAVE_DOXYGEN, test ! "x$DOXYGEN" = "x")
HAVE_DOT=NO
DOT_PATH=
AC_PATH_PROG(DOT, dot)
if ! test "x$DOT" = "x" ; then
   HAVE_DOT=YES
   DOT_PATH=`dirname $DOT`
fi
AC_SUBST([HAVE_DOT])
AC_SUBST([DOT_PATH])

dnl ------------------------------------------------------------------
dnl inatallation directory for libraries
dnl 
if test "x$prefix" != "xNONE" && test $ALICE_ROOT = $prefix ; then
libdir=\${exec_prefix}/lib/tgt_$ALICE_TARGET
AC_MSG_NOTICE([libdir set to $libdir])
fi
dnl ------------------------------------------------------------------
AC_CONFIG_FILES([Makefile])


AC_OUTPUT
dnl
dnl EOF
dnl

dnl AC_CONFIG_FILES([Makefile 
dnl 		 doc/Makefile
dnl 		 doc/doxygen.conf])
