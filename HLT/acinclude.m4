dnl
dnl $Id$
dnl
dnl  Copyright (C) 2002 Christian Holm Christensen <cholm@nbi.dk>
dnl
dnl  This library is free software; you can redistribute it and/or
dnl  modify it under the terms of the GNU Lesser General Public License
dnl  as published by the Free Software Foundation; either version 2.1
dnl  of the License, or (at your option) any later version.
dnl
dnl  This library is distributed in the hope that it will be useful,
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl  Lesser General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public
dnl  License along with this library; if not, write to the Free
dnl  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
dnl  02111-1307 USA
dnl
dnl ------------------------------------------------------------------
AC_DEFUN([AC_DEBUG],
[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_CXX])
  AC_MSG_CHECKING(whether to make debug objects)
  AC_ARG_ENABLE(debug,
    [AC_HELP_STRING([--enable-debug],[Enable debugging symbols in objects])],
    [],[enable_debug=no])
  if test "x$enable_debug" = "xno" ; then
    AC_DEFINE(NDEBUG)
    CFLAGS=`echo $CFLAGS | sed 's,-g,,'`
    CXXFLAGS=`echo $CXXFLAGS | sed 's,-g,,'`
  else
    AC_DEFINE(__DEBUG)
    AC_DEFINE(DEBUG)
    case $CXXFLAGS in
    *-g*) ;;
    *)    CXXFLAGS="$CXXFLAGS -g" ;;
    esac
    case $CFLAGS in
    *-g*) ;;
    *)    CFLAGS="$CFLAGS -g" ;;
    esac
  fi
  AC_MSG_RESULT($enable_debug 'CFLAGS=$CFLAGS')
])

dnl ------------------------------------------------------------------
AC_DEFUN([AC_OPTIMIZATION],
[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_CXX])

  AC_ARG_ENABLE(optimization,
    [AC_HELP_STRING([--disable-optimization],[Enable optimization of objects])],
    [],[enable_optimization=yes])

  AC_MSG_CHECKING(for optimiztion level)

  changequote(<<, >>)dnl
  if test "x$enable_optimization" = "xno" ; then
    CFLAGS=`echo   $CFLAGS   | sed 's,-O\([0-9][0-9]*\|\),,'`
    CXXFLAGS=`echo $CXXFLAGS | sed 's,-O\([0-9][0-9]*\|\),,'`
  elif test "x$enable_optimization" = "xyes" ; then
    case $CXXFLAGS in
    *-O*) ;;
    *)    CXXFLAGS="$CXXFLAGS -O2" ;;
    esac
    case $CFLAGS in
    *-O*) ;;
    *)    CFLAGS="$CXXFLAGS -O2" ;;
    esac
  else
    CFLAGS=`echo   $CFLAGS   | sed "s,-O\([0-9][0-9]*\|\),-O$enable_optimization,"`
    CXXFLAGS=`echo $CXXFLAGS | sed "s,-O\([0-9][0-9]*\|\),-O$enable_optimization,"`
  fi
  changequote([, ])dnl
  AC_MSG_RESULT($enable_optimization 'CFLAGS=$CFLAGS')
])

dnl ------------------------------------------------------------------

dnl
dnl Autoconf macro to check for existence or ROOT on the system
dnl Synopsis:
dnl
dnl  ROOT_PATH([MINIMUM-VERSION, [ACTION-IF-FOUND, [ACTION-IF-NOT-FOUND]]])
dnl
dnl Some examples: 
dnl 
dnl    ROOT_PATH(3.03/05, , AC_MSG_ERROR(Your ROOT version is too old))
dnl    ROOT_PATH(, AC_DEFINE([HAVE_ROOT]))
dnl 
dnl The macro defines the following substitution variables
dnl
dnl    ROOTCONF           full path to root-config
dnl    ROOTEXEC           full path to root
dnl    ROOTCINT           full path to rootcint
dnl    ROOTLIBDIR         Where the ROOT libraries are 
dnl    ROOTINCDIR         Where the ROOT headers are 
dnl    ROOTCFLAGS         Extra compiler flags
dnl    ROOTLIBS           ROOT basic libraries 
dnl    ROOTGLIBS          ROOT basic + GUI libraries
dnl    ROOTAUXLIBS        Auxilary libraries and linker flags for ROOT
dnl    ROOTAUXCFLAGS      Auxilary compiler flags 
dnl    ROOTRPATH          Same as ROOTLIBDIR
dnl
dnl The macro will fail if root-config and rootcint isn't found.
dnl
dnl Christian Holm Christensen <cholm@nbi.dk>
dnl
AC_DEFUN([ROOT_PATH],
[
  AC_ARG_WITH(rootsys,
  [  --with-rootsys          top of the ROOT installation directory],
    user_rootsys=$withval,
    user_rootsys="none")
  if test ! x"$user_rootsys" = xnone; then
    rootbin="$user_rootsys/bin"
  elif test ! x"$ROOTSYS" = x ; then 
    rootbin="$ROOTSYS/bin"
  else 
   rootbin=$PATH
  fi
  AC_PATH_PROG(ROOTCONF, root-config , no, $rootbin)
  AC_PATH_PROG(ROOTEXEC, root , no, $rootbin)
  AC_PATH_PROG(ROOTCINT, rootcint , no, $rootbin)
	
  if test ! x"$ROOTCONF" = "xno" && \
     test ! x"$ROOTCINT" = "xno" ; then 

    # define some variables 
    ROOTLIBDIR=`$ROOTCONF --libdir`
    ROOTINCDIR=`$ROOTCONF --incdir`
    ROOTCFLAGS=`$ROOTCONF --noauxcflags --cflags` 
    ROOTLIBS=`$ROOTCONF --noauxlibs --noldflags --libs`
    ROOTGLIBS=`$ROOTCONF --noauxlibs --noldflags --glibs`
    ROOTAUXCFLAGS=`$ROOTCONF --auxcflags`
    ROOTAUXLIBS=`$ROOTCONF --auxlibs`
    ROOTRPATH=$ROOTLIBDIR
	
    if test $1 ; then 
      AC_MSG_CHECKING(wether ROOT version >= [$1])
      vers=`$ROOTCONF --version | tr './' ' ' | awk 'BEGIN { FS = " "; } { printf "%d", ($''1 * 1000 + $''2) * 1000 + $''3;}'`
      requ=`echo $1 | tr './' ' ' | awk 'BEGIN { FS = " "; } { printf "%d", ($''1 * 1000 + $''2) * 1000 + $''3;}'`
      if test $vers -lt $requ ; then 
        AC_MSG_RESULT(no)
	no_root="yes"
      else 
        AC_MSG_RESULT(yes)
      fi
    fi
  else
    # otherwise, we say no_root
    no_root="yes"
  fi

  AC_SUBST(ROOTLIBDIR)
  AC_SUBST(ROOTINCDIR)
  AC_SUBST(ROOTCFLAGS)
  AC_SUBST(ROOTLIBS)
  AC_SUBST(ROOTGLIBS) 
  AC_SUBST(ROOTAUXLIBS)
  AC_SUBST(ROOTAUXCFLAGS)
  AC_SUBST(ROOTRPATH)

  if test "x$no_root" = "x" ; then 
    ifelse([$2], , :, [$2])     
  else 
    ifelse([$3], , :, [$3])     
  fi
])

dnl ------------------------------------------------------------------
dnl
dnl Autoconf macro to check conditions for an HLT module
dnl          - header dependencies
dnl          - library dependencies
dnl          - AliRoot availability
dnl The macro also exports the --enable/--disable option for the
dnl module.
dnl
dnl Synopsis:
dnl
dnl  ROOT_PATH([module],
dnl            [headers], [additional CPPFLAGS],
dnl            [libraries], [LD flags], [LIBS],
dnl            [circ libraries], [LD flags], [LIBS])
dnl
dnl First argument is the module name.
dnl
dnl Arg 2 and 3 allow to check a list of header files and to specify
dnl additional CPPFLAGS which might be necessary to perform the checks.
dnl 
dnl Arg 4 to 6 allow to check a list of libraries, with additional
dnl LDFLAGS and LIBS to be specified.
dnl
dnl Arg 7 to 9 is the same for libraries with circular dependencies
dnl among each other.
dnl
dnl Return:
dnl enable_module=yes
dnl   - AliRoot was found &&
dnl   - header files found &&
dnl   - module not disabled
dnl
dnl enable_module=no...requires.AliRoot
dnl   - AliRoot not found
dnl
dnl enable_module=missheader
dnl   - one of the specified header files is missing
dnl
dnl Libraries are probed and the variable ALIHLTMODULE_LIBS is set to
dnl all libraries which could be found.
dnl 
dnl Some examples:
dnl
dnl Matthias Richter <Matthias.Richter@uib.no>
AC_DEFUN([CHECK_HLTMODULE],
[
AH_TEMPLATE([HLT_[$1]],[hlt [$1] library])
AC_ARG_ENABLE([$1],
  [AC_HELP_STRING([--disable-[$1]],
      [   compile the $1 library ])],
  [if test "x$enableval" = "xno";
     then enable_module=no
   elif test "x$disable_all" = "xyes"; then
     # do checks if library has been enabled after global disable
     enable_module=yes
   else
     enable_module=force
  fi],
  [if test "x$disable_all" = "xyes"; then
   enable_module=no...modules.disabled
   else
   enable_module=yes
   fi # if test "x$disable_all" = "yes"
  ])
  if test "x$enable_module" = "xyes"; then
   if test "x$have_aliroot" = "xno" ; then
     enable_module="no...requires.AliRoot"
   else
     AC_MSG_NOTICE([-------------------------------------------------])
     AC_MSG_NOTICE([checking dependencies for [$1] library])

     AC_LANG_PUSH(C++)
     save_CPPFLAGS="$CPPFLAGS"
     save_LDFLAGS="$LDFLAGS"
     save_LIBS="$LIBS"
     CPPFLAGS="$save_CPPFLAGS [$3]"

     AC_CHECK_HEADERS([$2], [], [enable_module="missheader"])

     dnl ==========================================================================
     dnl
     dnl required header files and libraries for the AliHLTxxx library  
     dnl
     
     dnl ROOT/AliRoot libs needed by libAliHLTxxx
     CHECKLIBS="[$4]"
     ALIHLTMODULE_LIBS=
     LDFLAGS="$save_LDFLAGS [$5]"
     for CHECKLIB in $CHECKLIBS ; do
       LIBS="$save_LIBS $ROOTLIBS [$6] $ALIHLTMODULE_LIBS"
       AC_CHECK_LIB([$CHECKLIB],[_init], [ALIHLTMODULE_LIBS="$ALIHLTMODULE_LIBS -l$CHECKLIB"])
     done

     dnl libs with circular dependencies needed by libAliHLTxxx
     CHECKLIBS="[$7]"
     CIRCULARS=
     for dep in [$7]; do 
       CIRCULARS="$CIRCULARS -l$dep"
     done
     ALIHLTMODULE_LIBS="$ALIHLTMODULE_LIBS"
     LDFLAGS="$save_LDFLAGS [$8]"
     for CHECKLIB in $CHECKLIBS ; do
       CIRCULARS=`echo $CIRCULARS | sed -e "s|-l$CHECKLIB||"`
       LIBS="$save_LIBS [$9] $CIRCULARS $ALIHLTMODULE_LIBS"
       AC_CHECK_LIB([$CHECKLIB],[_init], [ALIHLTMODULE_LIBS="$ALIHLTMODULE_LIBS -l$CHECKLIB"])
     done
     CPPFLAGS="$save_CPPFLAGS"
     LDFLAGS="$save_LDFLAGS"
     LIBS="$save_LIBS"  
     AC_LANG_POP(C++)

   fi # if test "x$have_aliroot" = "xno"
  fi  # if test "x$enable_module" = "xyes"
])


#
# EOF
#
