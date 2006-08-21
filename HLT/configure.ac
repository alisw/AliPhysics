dnl -*- mode: autoconf -*- 
dnl
dnl $Id$
dnl template for the configuration script for the Alice HLT 
dnl framework and components
dnl 
dnl ------------------------------------------------------------------
AC_INIT([Alice High Level Trigger] , [0.2], [Matthias.Richter@ift.uib.no], alice-hlt)

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
CPPFLAGS="$CPPFLAGS -DROOTVERSION=\\\"`${ROOTCONF} --version`\\\""
CPPFLAGS="$CPPFLAGS -DALIROOTVERSION=\\\"Unknown\\\""
CPPFLAGS="$CPPFLAGS ${ROOTCFLAGS}"
LDFLAGS="$LDFLAGS -L${ROOTLIBDIR} ${ROOTLIBS}"

dnl ------------------------------------------------------------------
dnl check for AliRoot
dnl AC_LANG_PUSH(C++)
dnl AC_MSG_CHECKING([for AliRoot])
dnl have_aliroot=no
dnl ALIROOTLIBDIR=${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}
dnl ALIROOTINCDIR=${ALICE_ROOT}/include
dnl if test "x$ALICE_ROOT" != "x" \
dnl    && test -d ${ALIROOTLIBDIR} \
dnl    && test -d ${ALIROOTINCDIR}; then
dnl have_aliroot=yes
dnl else
dnl ALIROOTLIBDIR=
dnl ALIROOTINCDIR=
dnl fi
dnl AC_MSG_RESULT([$have_aliroot])

dnl if test ! "x$have_aliroot" = "x" ; then
dnl CPPFLAGS="$CPPFLAGS -I${ALIROOTINCDIR}"
dnl LDFLAGS="$LDFLAGS -L${ALIROOTLIBDIR} -lESD -lGeom -lMinuit -lVMC -lEG -lRAWData -ldl"
dnl have_alisteer=no
dnl AC_CHECK_LIB(STEER, [AliRun], [have_alisteer=yes])
dnl fi
dnl AC_LANG_POP(C++)


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
if test "x$enable_tpc" = "xyes" ; then 
  AC_DEFINE(HLT_TPC)
fi
AM_CONDITIONAL(EN_HLT_TPC, test x$enable_tpc = xyes)
AC_MSG_RESULT([$enable_tpc])

dnl ------------------------------------------------------------------
if test "x$enable_tpc" = "xyes" ; then 
  AC_MSG_CHECKING([whether to compile tpc reference library])
  AH_TEMPLATE([HLT_TPCREF],[hlt tpc reference library])
  AC_ARG_ENABLE(tpcref,
    [AC_HELP_STRING([--enable-tpcref],
        [compile the tpc reference library ])],
    [],[enable_tpcref=no])
  if test "x$enable_tpcref" = "xyes" ; then 
    AC_DEFINE(HLT_TPCREF)
  fi
  AM_CONDITIONAL(EN_HLT_TPCREF, test x$enable_tpcref = xyes)
  AC_MSG_RESULT([$enable_tpcref])
fi

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
AC_CONFIG_FILES([Makefile 
		 BASE/Makefile
		 SampleLib/Makefile
		 TPCLib/Makefile
		 TPCLib/Ref/Makefile
		 doc/Makefile
		 doc/doxygen.conf])


AC_OUTPUT
dnl
dnl EOF
dnl

