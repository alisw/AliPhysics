dnl -*- mode: autoconf -*- 
dnl
dnl $Id$
dnl template for the configuration script for the Alice HLT 
dnl framework and components
dnl 
dnl ------------------------------------------------------------------
AC_INIT([Alice High Level Trigger] , [0.5-dev], [Matthias.Richter@ift.uib.no], alice-hlt)

dnl ------------------------------------------------------------------
dnl the package from CVS contains the old Makefiles as well. In order to
dnl prevent them from becoming overwritten, we require a separate build
dnl directory
if test "`dirname $0`" = "." ; then
   AC_ERROR([please run the script from a separate build directory])
fi

dnl ------------------------------------------------------------------
AC_CANONICAL_SYSTEM
AC_PREFIX_DEFAULT(${PWD})
AC_CONFIG_SRCDIR(BASE/AliHLTComponent.cxx)
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_DEBUG
AC_OPTIMIZATION

dnl ------------------------------------------------------------------
dnl
dnl Check for ROOT
dnl
ROOT_PATH(, [have_root=1], [AC_ERROR([Stop! The HLT package needs ROOT.])])
AM_CONDITIONAL(HAVE_ROOT, test "x$have_root" = "x1" )

dnl ------------------------------------------------------------------
# TODO: make this configurable through arguments
#Define whether you want to run with ALIROOT or only ROOT
AH_TEMPLATE([ALIHLT_USEPACKAGE],[running environment])
ALIHLT_USEPACKAGE=ALIROOT
#ALIHLT_USEPACKAGE=ROOT
#ALIHLT_USEPACKAGE=STANDALONE
AC_DEFINE(use_aliroot)
AC_DEFINE(use_root)
CPPFLAGS="$CPPFLAGS ${ROOTCFLAGS}"

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
  # we check for the libSTEER library which is linked to
  # - the ROOT libs libGeom libMinuit libVMC libEG
  # - the AliRoot libESD libRAWData (libRAWDatarec from v4-04-Rev-07)
  # - ROOT libCint needs also libdl
  ALIROOT_CPPFLAGS="-I${ALIROOTINCDIR} -I${ALICE_ROOT}/RAW"
  ALIROOT_LDFLAGS="-L${ALIROOTLIBDIR}"
  ALIROOT_LIBS="-lESD -lGeom -lMinuit -lVMC -lEG -ldl"
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

  if test ! "x$have_aliroot" = "xno" ; then
    # the HLTTPCLib needs to link agains TPCbase and TPCrec
    # TPC library for AliTPCParam and AliSimDigits used in 
    # AliHLTTPCFileHandler.h and AliHLTTPCDisplay.cxx
    saveALIROOT_CPPFLAGS="$ALIROOT_CPPFLAGS"
    ALIROOT_CPPFLAGS="$saveALIROOT_CPPFLAGS -I${ALICE_ROOT}/TPC"
    CPPFLAGS="$save_CPPFLAGS $ALIROOT_CPPFLAGS"
    have_alitpc=yes
    AC_MSG_CHECKING(for required classes in TPC libraries)
    if test ! "x$have_alitpc" = "xno" ; then
      save_ALIROOT_LIBS=$ALIROOT_LIBS
      ALIROOT_LIBS="$save_ALIROOT_LIBS -lTPCbase -lTPCrec"
      LIBS="$save_LIBS $ROOTLIBS $ALIROOT_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <AliSimDigits.h>
  				       #include <AliTPCParam.h>
				       #include <AliTPCParamSR.h>
				       #include <AliTPCDigitsArray.h>
				       #include <AliTPCClustersArray.h>
				       #include <AliTPCcluster.h>     
				       #include <AliTPCClustersRow.h>], 
                                      [AliSimDigits dig;
  				       AliTPCParam param;
				       AliTPCParamSR paramsr;
				       AliTPCDigitsArray digarray;
				       AliTPCClustersArray clustarray;
				       AliTPCcluster clust;
				       AliTPCClustersRow row])],
                                      [AC_DEFINE(HAVE_ALITPC)
				       saveALIROOT_CPPFLAGS=$ALIROOT_CPPFLAGS
				       saveALIROOT_LIBS=$ALIROOT_LIBS], 
  				      [have_alitpc=no
				       ALIROOT_CPPFLAGS=$saveALIROOT_CPPFLAGS
				       ALIROOT_LIBS=$saveALIROOT_LIBS])
    fi
    AC_MSG_RESULT($have_alitpc)
  fi

  if test ! "x$have_alirawdata" = "xno" ; then
    # check for header files needed by the TPCLib/AliHLTTPCDigitReaderPacket
    # the header file has been moved to the TPC module
    CPPFLAGS="$save_CPPFLAGS $ALIROOT_CPPFLAGS"
    # AliRawReaderMemory.h AliRawReaderFile.h already checked
    #AC_CHECK_HEADERS([AliRawReaderMemory.h AliRawReaderFile.h AliTPCRawStream.h])
    AC_CHECK_HEADERS([AliTPCRawStream.h])
  fi

  have_tpc_mapping=no
  if test ! "x$have_aliroot" = "xno" ; then
    # the TPCLib/AliHLTTPCRawReader needs the actual Pad layot of the TPC 
    # which is stored in the TPC/mapping/Patchx.data files from
    # AliRoot version v4-04-Release (May 2006)
    AH_TEMPLATE([HAVE_TPC_MAPPING],[TPC pad mapping available])
    TPC_PAD_MAPPING_PATH=$have_aliroot/TPC/mapping
    AC_CHECK_FILE( [$TPC_PAD_MAPPING_PATH/Patch0.data], 
		   [have_tpc_mapping=yes
		    AC_SUBST(TPC_PAD_MAPPING_PATH)
		    AC_DEFINE(HAVE_TPC_MAPPING)],
		   [])
    
  fi
  AC_MSG_CHECKING([for TPC mapping layout])
  AM_CONDITIONAL(USE_TPC_MAPPING, test x$have_tpc_mapping = xyes)
  AC_MSG_RESULT([$have_tpc_mapping])
  
  if test "x$have_aliroot" = "xno" ; then
    ALIROOT_CPPFLAGS=
    ALIROOT_LDFLAGS=
    ALIROOT_LIBS=
  fi
  CPPFLAGS=$save_CPPFLAGS
  LDFLAGS=$save_LDFLAGS
  LIBS=$save_LIBS
  if test ! "x$have_aliroot" = "xno" ; then
    CPPFLAGS="$CPPFLAGS -DALIROOTVERSION=\\\"Unknown\\\""
  fi
fi
AC_LANG_POP(C++)
AC_SUBST([ALIROOT_CPPFLAGS])
AC_SUBST([ALIROOT_LDFLAGS])
AC_SUBST([ALIROOT_LIBS])
CPPFLAGS="$CPPFLAGS -DROOTVERSION=\\\"`${ROOTCONF} --version`\\\""

dnl ------------------------------------------------------------------
dnl check for the HLT PubSub Framework
AC_MSG_CHECKING([for HLT PubSub Framework])
AC_ARG_WITH(pubsub, [installation path of the HLT PubSub framework],
		    [],
		    [ test -n $ALIHLT_DC_DIR && with_pubsub=$ALIHLT_DC_DIR ])
if test -z $with_pubsub || ! test -d $with_pubsub ; then
  with_pubsub=no
fi
AC_MSG_RESULT([$with_pubsub])
HOMER_VERSION=1
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
    HOMER_CPPFLAGS="-I${with_pubsub}/include/HOMER -I${with_pubsub}/include/HOMER/reader"
  elif test -d ${with_pubsub}/include/Util/HOMER ; then
    # location for HLT Framework versions after Sep 2006
    HOMER_CPPFLAGS="-I${with_pubsub}/include/Util/HOMER"
  elif test -d ${with_pubsub}/src/Util/HOMER/include ; then
    # fall back if include files were not installed (versions after Sep 06)
    HOMER_CPPFLAGS="-I${with_pubsub}/src/Util/HOMER/include"
  else
    # fall back if include files were not installed (versions before Sep 06)
    HOMER_CPPFLAGS="-I${with_pubsub}/src/Util/HOMER/reader/include -I${with_pubsub}/src/Util/HOMER/data/include"
  fi
  HOMER_LDFLAGS="-L${with_pubsub}/lib/`uname -s`-`uname -m`"
  CPPFLAGS="$save_CPPFLAGS $HOMER_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS $HOMER_LDFLAGS"
  LIBS="$save_LIBS "
  with_homer=no
  AC_LANG_PUSH(C++)
  AC_CHECK_HEADER([HOMERReader.h],
    [ AC_CHECK_LIB([HOMER],[_init],
      [with_homer=yes
       HOMER_LIBS="-lHOMER"],
      [ AC_CHECK_LIB([HOMERReader],[_init],
                     [with_homer=yes
		      HOMER_LIBS="-lHOMERReader"])])
      LIBS="$LIBS $HOMER_LIBS"
      AC_MSG_CHECKING([for HOMER version])
      dnl The Homer library has no versioning, so we do our own
      dnl version           description
      dnl ----------------------------------------------------------------------
      dnl   1          inintial version
      dnl   2          GetBlockStatusFlags added to the interface
      AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <HOMERReader.h>],
                                      [HOMERReader reader("dummy", 0);
 				       reader.GetBlockStatusFlags(0);])],
                                      [if test $HOMER_VERSION -lt 2 ; then HOMER_VERSION=2; fi],
 				      [])
      AC_MSG_RESULT([$HOMER_VERSION])
  ]) #AC_CHECK_HEADER([HOMERReader.h])
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
fi
AC_SUBST([HOMER_CPPFLAGS])
AC_SUBST([HOMER_LDFLAGS])
AC_SUBST([HOMER_LIBS])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to impose strict coding conventions])
AC_ARG_ENABLE(strict,
  [AC_HELP_STRING([--enable-strict],
      [enable coding convention checks ])],
  [],[enable_strict=no])
if test "x$enable_strict" = "xyes" ; then
   CPPFLAGS="$CPPFLAGS -Weffc++"
fi
AC_MSG_RESULT([$enable_strict])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to compile sample library])
AH_TEMPLATE([HLT_SAMPLE],[hlt sample library])
AC_ARG_ENABLE(sample,
  [AC_HELP_STRING([--enable-sample],
      [compile the sample library ])],
  [],[enable_sample=yes])
if test "x$enable_sample" = "xyes" ; then 
  AC_DEFINE(HLT_SAMPLE)
fi
AM_CONDITIONAL(EN_HLT_SAMPLE, test x$enable_sample = xyes)
AC_MSG_RESULT([$enable_sample])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to compile tpc library])
AH_TEMPLATE([HLT_TPC],[hlt tpc library])
AC_ARG_ENABLE(tpc,
  [AC_HELP_STRING([--enable-tpc],
      [compile the tpc library ])],
  [],[enable_tpc=yes])
if test "x$have_aliroot" = "xno" ; then
   enable_tpc="no...requires.AliRoot"
elif test "x$have_alitpc" = "xno" ; then
   enable_tpc="no...requires.AliRoot.TPC.libraries"
fi
if test "x$enable_tpc" = "xyes" ; then 
  AC_DEFINE(HLT_TPC)
fi
AM_CONDITIONAL(EN_HLT_TPC, test x$enable_tpc = xyes)
AC_MSG_RESULT([$enable_tpc])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to compile phos library])
AH_TEMPLATE([HLT_PHOS],[hlt phos library])
AC_ARG_ENABLE(phos,
  [AC_HELP_STRING([--enable-phos],
      [compile the phos library ])],
  [],[enable_phos=yes])
if test "x$have_aliroot" = "xno" ; then
   enable_phos="no...requires.AliRoot"
fi
if test "x$enable_phos" = "xyes" ; then 
  AC_DEFINE(HLT_PHOS)
fi
AM_CONDITIONAL(EN_HLT_PHOS, test x$enable_phos = xyes)
AC_MSG_RESULT([$enable_phos])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to enable HLT logging])
AH_TEMPLATE([NOLOGGING],[disable hlt logging])
AC_ARG_ENABLE(logging,
  [AC_HELP_STRING([--enable-logging],
      [enable logging])],
  [],[enable_logging=yes])
if test "x$enable_logging" != "xyes" ; then 
  AC_DEFINE(NOLOGGING)
fi
AC_MSG_RESULT([$enable_logging])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to enable saving MC data through the chain])
AH_TEMPLATE([DOMC],[MC saving])
AC_ARG_ENABLE(mc-saving,
  [AC_HELP_STRING([--enable-mc-saving],
      [enable saving MC data through the chain])],
  [],[enable_mc_saving=no])
if test "x$enable_mc_saving" = "xyes" ; then 
  AC_DEFINE(DOMC)
fi
AC_MSG_RESULT([$enable_mc_saving])

dnl ------------------------------------------------------------------
AC_MSG_CHECKING([whether to enable AliRoot NEWIO])
AH_TEMPLATE([use_newio],[AliRoot NEWIO])
AC_ARG_ENABLE(newio,
  [AC_HELP_STRING([--enable-newio],
      [enable AliRoot NEWIO ])],
  [],[enable_newio=yes])
if test "x$enable_newio" = "xyes" ; then 
  AC_DEFINE(use_newio)
fi
AC_MSG_RESULT([$enable_newio])

dnl ------------------------------------------------------------------
dnl certainly something old, but we keep the define
AC_MSG_CHECKING([whether to use ROWHOUGH])
AH_TEMPLATE([USEROWHOUGH],[HLT ROWHOUGH])
AC_ARG_ENABLE(rowhough,
  [AC_HELP_STRING([--enable-rowhough],
      [use ROWHOUGH ])],
  [],[enable_rowhough=no])
if test "x$enable_rowhough" = "xyes" ; then 
  AC_DEFINE(USEROWHOUGH)
fi
AC_MSG_RESULT([$enable_rowhough])

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
AC_CONFIG_FILES([Makefile 
		 BASE/Makefile
		 SampleLib/Makefile
		 TPCLib/Makefile
		 TPCLib/mapping2array.cxx
		 TPCLib/OnlineDisplay/Makefile
		 PHOS/Makefile
		 doc/Makefile
		 doc/doxygen.conf])


AC_OUTPUT
dnl
dnl EOF
dnl

