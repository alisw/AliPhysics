dnl -*- mode: autoconf -*- 
dnl
dnl $Id$
dnl template for the configuration script for the Alice HLT 
dnl Online Monitoring GUI
dnl 
dnl ------------------------------------------------------------------
AC_INIT([Alice HLT Online Monitoring GUI] , [1.5.3], [thaeder@kip.uni-heidelberg.de], AliHLTGUI)

dnl ------------------------------------------------------------------
AC_CANONICAL_SYSTEM
AC_PREFIX_DEFAULT(${PWD})
AC_CONFIG_SRCDIR(src/AliHLTGUI.cxx)
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_DEBUG
AC_OPTIMIZATION

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to apply strict checking])
AC_ARG_ENABLE(strict,
  [AC_HELP_STRING([--disable-strict],
      [disable the strict checking of dependencies during configure])],
  [],[enable_strict=yes])
AC_MSG_RESULT([$enable_strict])

dnl ------------------------------------------------------------------
dnl
dnl Check for Trolltechs Qt
dnl
BNV_HAVE_QT
test "x$have_qt" = "xno" && AC_ERROR([Stop! The package needs Qt])
AC_SUBST([QT_LIBS])
AC_SUBST([QT_DIR])

dnl ------------------------------------------------------------------
dnl
dnl Check for ROOT
dnl
ROOT_PATH(5.08, [have_root=1], [AC_ERROR([Stop! The package needs ROOT 5.08 or later.])])
AM_CONDITIONAL(HAVE_ROOT, test "x$have_root" = "x1" )
AC_SUBST([ROOTSYS])
ROOTBINDIR=`dirname $ROOTEXEC`
AC_SUBST([ROOTBINDIR])

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
    CPPFLAGS="$save_CPPFLAGS ${ROOTCFLAGS}"
    LDFLAGS="$save_LDFLAGS -L${ROOTLIBDIR}"
    LIBS="$save_LIBS $ROOTLIBS $ADD_ROOTLIBS"
    AC_CHECK_LIB([$CHECKLIB],[_init], [ADD_ROOTLIBS="$ADD_ROOTLIBS -l$CHECKLIB"])
  done

dnl in order to have the configure.ac somewhat similar to the AliRoot HLT module's
dnl one we keep this three lines commented. The build system of this package
dnl might be changed later in order to pu the ROOT CPPFLAGS, LDFLAGS and LIBS
dnl directly into the Makefile.am
dnl   CPPFLAGS=$save_CPPFLAGS
dnl   LDFLAGS=$save_LDFLAGS
dnl   LIBS=$save_LIBS
fi


dnl ------------------------------------------------------------------
dnl check for AliRoot features
AC_LANG_PUSH(C++)
have_aliroot=no
AC_ARG_WITH(aliroot,[  --with-aliroot	top of the AliRoot installation],
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
  ALIROOTBINDIR=${ALICE_ROOT}/bin/tgt_${ALICE_TARGET}
  ALIROOTLIBDIR=${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}
  ALIROOTINCDIR=${ALICE_ROOT}/include
  test -d ${ALIROOTBINDIR} || AC_MSG_WARN([can not find AliRoot binary directory $ALIROOTBINDIR])
  test -d ${ALIROOTLIBDIR} || AC_MSG_WARN([can not find AliRoot library directory $ALIROOTLIBDIR])
  test -d ${ALIROOTINCDIR} || AC_MSG_WARN([can not find AliRoot include directory $ALIROOTINCDIR])
fi

AC_MSG_CHECKING([for AliRoot])
if test "x$ALICE_ROOT" != "x" \
   && test -d ${ALIROOTBINDIR} \
   && test -d ${ALIROOTLIBDIR} \
   && test -d ${ALIROOTINCDIR}; then
  have_aliroot=$ALICE_ROOT
else
  ALIROOTBINDIR=
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
  # libSTEERBase present since Aug 7 2007
  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  CHECKLIB=STEERBase
  AC_CHECK_LIB([$CHECKLIB],[_init],[ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"])

  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  CHECKLIB=AOD
  AC_CHECK_LIB([$CHECKLIB],[_init],[ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"])

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
 				      [have_alirawdata=no
				       have_aliroot=no])
	]) dnl AC_CHECK_LIB RAWDatabase
  AC_MSG_CHECKING(for AliRawReader classes in RAWData libraries)
  AC_MSG_RESULT($have_alirawdata)

  LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
  CHECKLIB=STEER
  AC_CHECK_LIB([$CHECKLIB],[_init],[ALIROOT_LIBS="$ALIROOT_LIBS -l$CHECKLIB"],
				   [have_aliroot=no])

  if test "x$have_aliroot" = "xno" ; then
    ALIROOT_CPPFLAGS=
    ALIROOT_LDFLAGS=
    ALIROOT_LIBS=
  fi
  CPPFLAGS=$save_CPPFLAGS
  LDFLAGS=$save_LDFLAGS
  LIBS=$save_LIBS
fi

test "x$have_aliroot" = "xno" && \
     test "x$enable_strict" != "xno" && \
     AC_ERROR([Stop! The package needs AliRoot libraries. Check your setup])

AC_SUBST([ALIROOTLIBDIR])
AC_SUBST([ALIROOT_CPPFLAGS])
AC_SUBST([ALIROOT_LDFLAGS])
AC_SUBST([ALIROOTBINDIR])
AC_SUBST([ALIROOTLIBDIR])
AC_SUBST([ALIROOT_LIBS])

dnl ------------------------------------------------------------------
dnl Check for the AliHLTTPCLibrary
dnl
AC_ARG_WITH(alihlt, [  --with-alihlt		top of the AliRoot HLT module],
		    [],
		    [ test -n $ALIHLT_PREFIX && with_alihlt=$ALIHLT_PREFIX])
if ! test -z ${with_alihlt} \
  && test -d ${with_alihlt} \
  && test "x${with_alihlt}" != "x$ALICE_ROOT" ; then
  ALIHLT_INCDIRS="${with_alihlt}/include/alice-hlt"
  ALIHLT_LIBDIR=${with_alihlt}/lib
else
  ALIHLT_INCDIRS="${ALICE_ROOT}/HLT/BASE ${ALICE_ROOT}/HLT/BASE/HOMER"
  ALIHLT_LIBDIR=${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}
fi
ALIHLT_CPPFLAGS=`for i in ${ALIHLT_INCDIRS}; do echo -I$i ; done`
ALIHLT_LDFLAGS="-L$ALIHLT_LIBDIR"

ALIHLT_LIBS='-lHLTbase -lAliHLTTPC -lAliHLTTPCDisplay'
save_CPPFLAGS=$CPPFLAGS
save_LDFLAGS=$LDFLAGS
save_LIBS=$LIBS
CPPFLAGS="$ALIHLT_CPPFLAGS $save_CPPFLAGS $ALIROOT_CPPFLAGS"
LDFLAGS="$ALIHLT_LDFLAGS $save_LDFLAGS $ALIROOT_LDFLAGS $HOMER_LDFLAGS"
LIBS="$save_LIBS $ALIHLT_LIBS $ALIROOT_LIBS $HOMER_LIBS"
export LD_LIBRARY_PATH="$ALIHLT_LIBDIR:$HOMER_LIBDIR:$LD_LIBRARY_PATH"

have_alitpcdisplay=yes
ALIHLTTPCDISPLAY_VERSION=0
AC_CHECK_LIB([AliHLTTPCDisplay],[__libAliHLTTPCDisplayProperties],
  [
#    AC_MSG_NOTICE([$LD_LIBRARY_PATH])
    AC_MSG_CHECKING(AliHLTTPCDisplay library for HOMER support)
    AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <string.h>
				     #include <stdio.h>
				     extern const char* __libAliHLTTPCDisplayProperties;],
                                    [FILE* fp=fopen("conftest.out", "w");
				     if (fp) {fprintf(fp, __libAliHLTTPCDisplayProperties);fclose(fp);}
				     // return 0 (success) if string was found (!=NULL)
				     return (strstr(__libAliHLTTPCDisplayProperties, "HOMER")==NULL)])],
                                    [test -f conftest.out && echo -n "library properties: `cat conftest.out` "
				     ALIHLTTPCDISPLAY_VERSION=1],
				    [have_alitpcdisplay=no
				     test -f conftest.out && echo "library properties: `cat conftest.out`"
				     test "x$enable_strict" != "xno" && AC_MSG_ERROR([Stop! The AliHLTTPCDisplay library must be compiled with HOMER support])])
  ],
  [
    AC_MSG_CHECKING(for usability of the AliHLTTPCDisplay library)
    if test ! "x$have_alitpcdisplay" = "xno" ; then
      AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <AliHLTTPCDisplayMain.h>], 
                                      [AliHLTTPCDisplayMain display(NULL,NULL)])],
                                      [], 
				      [have_alitpcdisplay=no])
    fi
  ])
AC_MSG_RESULT($have_alitpcdisplay)
test "x$have_alitpcdisplay" != "xyes" && \
     test "x$enable_strict" != "xno" && \
     AC_ERROR([Stop! The package needs the AliHLTTPCDisplay library.
     Use --with-alihlt option to point to the AliRoot HLT module installation])

AC_MSG_CHECKING([version of AliHLTTPCDisplay library])
if test $ALIHLTTPCDISPLAY_VERSION -gt 0 ; then
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <AliHLTTPCDisplayMain.h>], 
                                     [AliHLTTPCDisplayMain dm(NULL,NULL) ;
 				       dm.SetZeroSuppression(0);
				       dm.GetNTimeBins();
				       dm.SetTimeBinMinMax(0,1);
				       dm.SetFrontDataSwitch(0);
				       dm.Set3DRawSwitch(0)])],
                                     [ALIHLTTPCDISPLAY_VERSION=2],
 				     [])
fi
AC_MSG_RESULT([$ALIHLTTPCDISPLAY_VERSION])

CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS
LIBS=$save_LIBS

AC_SUBST([ALIHLT_INCDIRS])
AC_SUBST([ALIHLT_LIBDIR])
AC_SUBST([ALIHLT_CPPFLAGS])
AC_SUBST([ALIHLT_LDFLAGS])
AC_SUBST([ALIHLT_LIBS])
AC_SUBST([ALIHLTTPCDISPLAY_VERSION])

AC_LANG_POP(C++)

dnl ------------------------------------------------------------------
dnl check for the HLT PubSub Framework
dnl namely for the existence of the HOMER library
dnl from Sep 2007, the HOMER lib has been incorporated into the alice-hlt
dnl package. It os though possible to choose an external. library
dnl In order to make the origni of the HOMER lib clear, the one in AliRoot
dnl got the name libAliHLTHOMER
AC_MSG_CHECKING([for HLT PubSub Framework])
AC_ARG_WITH(pubsub, [installation path of the HLT PubSub framework],
		    [],
		    [ test -n $ALIHLT_DC_DIR && with_pubsub=$ALIHLT_DC_DIR ])
if test -z $with_pubsub || ! test -d $with_pubsub ; then
  with_pubsub=no
fi
AC_MSG_RESULT([$with_pubsub])
HOMER_VERSION=2
HOMER_LIBS=

AH_TEMPLATE([HAVE_HOMERREADER],[the HLT PubSub Homer Reader interface])
if test "x$with_pubsub" != "xno" ; then
  save_CPPFLAGS=$CPPFLAGS
  save_LDFLAGS=$LDFLAGS
  save_LIBS=$LIBS
  # currently the different versions of the HLT PubSub framework have a different
  # directory layout
  if test -d ${with_pubsub}/include/HOMER ; then
    # the 'early' location of the include files with separated HOMER sub dirs
    HOMER_INCDIRS="${with_pubsub}/include/HOMER ${with_pubsub}/include/HOMER/reader"
  elif test -d ${with_pubsub}/include/Util/HOMER ; then
    # location for HLT Framework versions after Sep 2006
    HOMER_INCDIRS="${with_pubsub}/include/Util/HOMER"
  elif test -d ${with_pubsub}/src/Util/HOMER/include ; then
    # fall back if include files were not installed (versions after Sep 06)
    HOMER_INCDIRS="${with_pubsub}/src/Util/HOMER/include"
  else
    # fall back if include files were not installed (versions before Sep 06)
    HOMER_INCDIRS="${with_pubsub}/src/Util/HOMER/reader/include ${with_pubsub}/src/Util/HOMER/data/include"
  fi
  HOMER_CPPFLAGS=`for i in ${HOMER_INCDIRS}; do echo -n "-I$i " ; done`

  AC_ARG_ENABLE(pubsub-debug,
    [AC_HELP_STRING([--disable-pubsub-debug],
        [force the production version of the PubSub framework])],
    [],[enable_pubsub_debug=yes])
  homer_type="none"
  HOMER_TARGET="`uname -s`-`uname -m`"
  AC_MSG_CHECKING([for HLT PubSub Framework release type])
  if test -d "${with_pubsub}/lib/${HOMER_TARGET}-release" ; then
    homer_type="production"
    HOMER_TARGET="${HOMER_TARGET}-release"
  elif test "x$enable_pubsub_debug" = "xyes" && test -d "${with_pubsub}/lib/${HOMER_TARGET}-debug" ; then
    homer_type="debug"
    HOMER_TARGET="${HOMER_TARGET}-debug"
  fi
  AC_MSG_RESULT([${homer_type}])
  HOMERREADER_HEADER=HOMERReader.h
  HOMER_BINDIR="${with_pubsub}/bin/${HOMER_TARGET}"
  HOMER_LIBDIR="${with_pubsub}/lib/${HOMER_TARGET}"
  HOMER_LDFLAGS="-L${HOMER_LIBDIR}"
  CHECK_LIB=HOMER
  CPPFLAGS="$save_CPPFLAGS $HOMER_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS $HOMER_LDFLAGS"
  LIBS="$save_LIBS "
else
  HOMERREADER_HEADER=AliHLTHOMERReader.h
  HOMER_BINDIR=
  HOMER_LIBDIR=$ALIHLT_LIBDIR
  HOMER_LDFLAGS="-L${HOMER_LIBDIR}"
  CHECK_LIB=AliHLTHOMER
  HOMER_CPPFLAGS=`for i in ${ALIHLT_INCDIRS}; do echo -I$i ; done`
  CPPFLAGS="$save_CPPFLAGS $HOMER_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS $HOMER_LDFLAGS"
  LIBS="$save_LIBS "
fi
  with_homer=no
  AC_LANG_PUSH(C++)
  AC_CHECK_HEADER([$HOMERREADER_HEADER],
    [ AC_CHECK_LIB([$CHECK_LIB],[_init],
      [with_homer=yes
       HOMER_LIBS="-l$CHECK_LIB"],
      [ AC_CHECK_LIB([HOMERReader],[_init],
                     [with_homer=yes
		      HOMER_LIBS="-lHOMERReader"])])
      LIBS="$LIBS $HOMER_LIBS"
      AC_MSG_CHECKING([version of HOMER library])
      dnl The Homer library has no versioning, so we do our own
      dnl version           description
      dnl ----------------------------------------------------------------------
      dnl   1          inintial version
      dnl   2          GetBlockStatusFlags added to the interface
      AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <$HOMERREADER_HEADER>],
                                      [HOMERReader reader((const void*)0, (int)0);
 				       reader.GetBlockStatusFlags(0);])],
                                      [],
 				      [HOMER_VERSION=1])
      AC_MSG_RESULT([$HOMER_VERSION])
  ],
  [HOMERREADER_HEADER=]) #AC_CHECK_HEADER([$HOMERREADER_HEADER])
  AC_LANG_POP(C++)
  if test "x$with_homer" = "xno" ; then
    HOMER_CPPFLAGS=
    HOMER_LDFLAGS=
    HOMER_LIBS=
  else
    HOMER_CPPFLAGS="$HOMER_CPPFLAGS -DHOMER_VERSION=$HOMER_VERSION"
    AC_DEFINE(HAVE_HOMERREADER)
  fi
  CPPFLAGS="$save_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS"
  LIBS="$save_LIBS"  

test "x$with_homer" != "xyes" && \
     test "x$enable_strict" != "xno" && \
     AC_ERROR([Stop! The package needs the HOMER library.
      Use  --with-pubsub option to point to the HLT PubSub installation])

AC_SUBST([HOMER_CPPFLAGS])
AC_SUBST([HOMER_LDFLAGS])
AC_SUBST([HOMER_LIBDIR])
AC_SUBST([HOMER_BINDIR])
AC_SUBST([HOMER_INCDIRS])
AC_SUBST([HOMER_LIBS])


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
AC_CONFIG_FILES([Makefile
		 doc/Makefile
		 doc/doxygen.conf
		 src/startAliHLTGUI
		 src/flags])


AC_OUTPUT
dnl
dnl EOF
dnl

