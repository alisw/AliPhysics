dnl
dnl $Id$
dnl
dnl This files contains a couple of autoconf macros, see macro 
dnl description for copyright notice.
dnl
dnl
dnl ------------------------------------------------------------------
dnl
dnl AC_DEBUG and AC_OPTIMIZATION
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
    [],[enable_debug=yes])
  if test "x$enable_debug" = "xno" ; then
    CFLAGS=`echo $CFLAGS | sed 's,-g,,'`
    CXXFLAGS=`echo $CXXFLAGS | sed 's,-g,,'`
  else
    AC_DEFINE(__DEBUG)
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
    [AC_HELP_STRING([--enable-optimization],[Enable optimization of objects])],
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
dnl ROOT_PATH
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
dnl BNV_HAVE_QT taken from the autoconf macro archive
dnl http://www.gnu.org/software/ac-archive/
dnl
dnl @synopsis BNV_HAVE_QT [--with-Qt-dir=DIR] [--with-Qt-lib-dir=DIR] [--with-Qt-lib=LIB]
dnl @synopsis BNV_HAVE_QT [--with-Qt-include-dir=DIR] [--with-Qt-bin-dir=DIR] [--with-Qt-lib-dir=DIR] [--with-Qt-lib=LIB]
dnl
dnl @summary Search for Trolltech's Qt GUI framework.
dnl
dnl Searches common directories for Qt include files, libraries and Qt
dnl binary utilities. The macro supports several different versions of
dnl the Qt framework being installed on the same machine. Without
dnl options, the macro is designed to look for the latest library,
dnl i.e., the highest definition of QT_VERSION in qglobal.h. By use of
dnl one or more options a different library may be selected. There are
dnl two different sets of options. Both sets contain the option
dnl --with-Qt-lib=LIB which can be used to force the use of a
dnl particular version of the library file when more than one are
dnl available. LIB must be in the form as it would appear behind the
dnl "-l" option to the compiler. Examples for LIB would be "qt-mt" for
dnl the multi-threaded version and "qt" for the regular version. In
dnl addition to this, the first set consists of an option
dnl --with-Qt-dir=DIR which can be used when the installation conforms
dnl to Trolltech's standard installation, which means that header files
dnl are in DIR/include, binary utilities are in DIR/bin and the library
dnl is in DIR/lib. The second set of options can be used to indicate
dnl individual locations for the header files, the binary utilities and
dnl the library file, in addition to the specific version of the
dnl library file.
dnl
dnl The following shell variable is set to either "yes" or "no":
dnl
dnl   have_qt
dnl
dnl Additionally, the following variables are exported:
dnl
dnl   QT_CXXFLAGS
dnl   QT_LIBS
dnl   QT_MOC
dnl   QT_UIC
dnl   QT_DIR
dnl
dnl which respectively contain an "-I" flag pointing to the Qt include
dnl directory (and "-DQT_THREAD_SUPPORT" when LIB is "qt-mt"), link
dnl flags necessary to link with Qt and X, the name of the meta object
dnl compiler and the user interface compiler both with full path, and
dnl finaly the variable QTDIR as Trolltech likes to see it defined (if
dnl possible).
dnl
dnl Example lines for Makefile.in:
dnl
dnl   CXXFLAGS = @QT_CXXFLAGS@
dnl   MOC      = @QT_MOC@
dnl
dnl After the variables have been set, a trial compile and link is
dnl performed to check the correct functioning of the meta object
dnl compiler. This test may fail when the different detected elements
dnl stem from different releases of the Qt framework. In that case, an
dnl error message is emitted and configure stops.
dnl
dnl No common variables such as $LIBS or $CFLAGS are polluted.
dnl
dnl Options:
dnl
dnl --with-Qt-dir=DIR: DIR is equal to $QTDIR if you have followed the
dnl installation instructions of Trolltech. Header files are in
dnl DIR/include, binary utilities are in DIR/bin and the library is in
dnl DIR/lib.
dnl
dnl --with-Qt-include-dir=DIR: Qt header files are in DIR.
dnl
dnl --with-Qt-bin-dir=DIR: Qt utilities such as moc and uic are in DIR.
dnl
dnl --with-Qt-lib-dir=DIR: The Qt library is in DIR.
dnl
dnl --with-Qt-lib=LIB: Use -lLIB to link with the Qt library.
dnl
dnl If some option "=no" or, equivalently, a --without-Qt-* version is
dnl given in stead of a --with-Qt-*, "have_qt" is set to "no" and the
dnl other variables are set to the empty string.
dnl
dnl @category InstalledPackages
dnl @author Bastiaan Veelo <Bastiaan@Veelo.net>
dnl @version 2006-03-12
dnl @license AllPermissive

dnl Copyright (C) 2001, 2002, 2003, 2005, 2006 Bastiaan Veelo

dnl THANKS! This code includes bug fixes and contributions made by:
dnl Tim McClarren,
dnl Dennis R. Weilert,
dnl Qingning Huo,
dnl Brian Mingus,
dnl Jens Hannemann,
dnl Pavel Roskin,
dnl Scott J. Bertin.

dnl ChangeLog
dnl 2006-03-12  * Hide output of ls and fix an m4 quoting problem (due to Scott J. Bertin).
dnl 2006-02-13  * Check compiler return value instead of parsing the error stream,
dnl               which detected warnings as false negatives (due to Jens Hannemann).
dnl 2006-02-02  * Spelling of "Success".
dnl             * Fixed unsave test for $bnv_qt_lib without quotes.
dnl             * Put dnl in front of all comments.
dnl             * Changed -l$bnv_qt_lib_dir into -L$bnv_qt_lib_dir (all due to Pavel Roskin).
dnl 2006-01-19  * Support for 64bit architectures.
dnl             * Updated documentation.
dnl 2006-01-18: * Fix "cat: bnv_qt_test.c: No such file or directory" (due to Jens Hannemann).
dnl             * Hide output of failing ls.
dnl 2006-01-11: * Check in /Developer on Mac OS X; Check in $QTDIR (due to Brian Mingus).

dnl Calls BNV_PATH_QT_DIRECT (contained in this file) as a subroutine.
AC_DEFUN([BNV_HAVE_QT],
[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([AC_PATH_X])
  AC_REQUIRE([AC_PATH_XTRA])

  AC_MSG_CHECKING(for Qt)

  AC_ARG_WITH([Qt-dir],
    [  --with-Qt-dir=DIR       DIR is equal to $QTDIR if you have followed the
                          installation instructions of Trolltech. Header
                          files are in DIR/include, binary utilities are
                          in DIR/bin. The library is in DIR/lib, unless
			  --with-Qt-lib-dir is also set.])
  AC_ARG_WITH([Qt-include-dir],
    [  --with-Qt-include-dir=DIR
                          Qt header files are in DIR])
  AC_ARG_WITH([Qt-bin-dir],
    [  --with-Qt-bin-dir=DIR   Qt utilities such as moc and uic are in DIR])
  AC_ARG_WITH([Qt-lib-dir],
    [  --with-Qt-lib-dir=DIR   The Qt library is in DIR])
  AC_ARG_WITH([Qt-lib],
    [  --with-Qt-lib=LIB       Use -lLIB to link with the Qt library])
  if test x"$with_Qt_dir" = x"no" ||
     test x"$with_Qt_include-dir" = x"no" ||
     test x"$with_Qt_bin_dir" = x"no" ||
     test x"$with_Qt_lib_dir" = x"no" ||
     test x"$with_Qt_lib" = x"no"; then
    # user disabled Qt. Leave cache alone.
    have_qt="User disabled Qt."
  else
    # "yes" is a bogus option
    if test x"$with_Qt_dir" = xyes; then
      with_Qt_dir=
    fi
    if test x"$with_Qt_include_dir" = xyes; then
      with_Qt_include_dir=
    fi
    if test x"$with_Qt_bin_dir" = xyes; then
      with_Qt_bin_dir=
    fi
    if test x"$with_Qt_lib_dir" = xyes; then
      with_Qt_lib_dir=
    fi
    if test x"$with_Qt_lib" = xyes; then
      with_Qt_lib=
    fi
    # No Qt unless we discover otherwise
    have_qt=no
    # Check whether we are requested to link with a specific version
    if test x"$with_Qt_lib" != x; then
      bnv_qt_lib="$with_Qt_lib"
    fi
    # Check whether we were supplied with an answer already
    if test x"$with_Qt_dir" != x; then
      have_qt=yes
      bnv_qt_dir="$with_Qt_dir"
      bnv_qt_include_dir="$with_Qt_dir/include"
      bnv_qt_bin_dir="$with_Qt_dir/bin"
      bnv_qt_lib_dir="$with_Qt_dir/lib"
      # Only search for the lib if the user did not define one already
      if test x"$bnv_qt_lib" = x; then
        bnv_qt_lib="`ls $bnv_qt_lib_dir/libqt* | sed -n 1p |
                     sed s@$bnv_qt_lib_dir/lib@@ | [sed s@[.].*@@]`"
      fi
      bnv_qt_LIBS="-L$bnv_qt_lib_dir -l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
    else
      # Use cached value or do search, starting with suggestions from
      # the command line
      AC_CACHE_VAL(bnv_cv_have_qt,
      [
        # We are not given a solution and there is no cached value.
        bnv_qt_dir=NO
        bnv_qt_include_dir=NO
        bnv_qt_lib_dir=NO
        if test x"$bnv_qt_lib" = x; then
          bnv_qt_lib=NO
        fi
        BNV_PATH_QT_DIRECT
        if test "$bnv_qt_dir" = NO ||
           test "$bnv_qt_include_dir" = NO ||
           test "$bnv_qt_lib_dir" = NO ||
           test "$bnv_qt_lib" = NO; then
          # Problem with finding complete Qt.  Cache the known absence of Qt.
          bnv_cv_have_qt="have_qt=no"
        else
          # Record where we found Qt for the cache.
          bnv_cv_have_qt="have_qt=yes                  \
                       bnv_qt_dir=$bnv_qt_dir          \
               bnv_qt_include_dir=$bnv_qt_include_dir  \
                   bnv_qt_bin_dir=$bnv_qt_bin_dir      \
                      bnv_qt_LIBS=\"$bnv_qt_LIBS\""
        fi
      ])dnl
      eval "$bnv_cv_have_qt"
    fi # all $bnv_qt_* are set
  fi   # $have_qt reflects the system status
  if test x"$have_qt" = xyes; then
    QT_CXXFLAGS="-I$bnv_qt_include_dir"
    if test x"$bnv_qt_lib" = xqt-mt; then
        QT_CXXFLAGS="$QT_CXXFLAGS -DQT_THREAD_SUPPORT"
    fi
    QT_DIR="$bnv_qt_dir"
    QT_LIBS="$bnv_qt_LIBS"
    # If bnv_qt_dir is defined, utilities are expected to be in the
    # bin subdirectory
    if test x"$bnv_qt_dir" != x; then
        if test -x "$bnv_qt_dir/bin/uic"; then
          QT_UIC="$bnv_qt_dir/bin/uic"
        else
          # Old versions of Qt don't have uic
          QT_UIC=
        fi
      QT_MOC="$bnv_qt_dir/bin/moc"
    else
      # Or maybe we are told where to look for the utilities
      if test x"$bnv_qt_bin_dir" != x; then
        if test -x "$bnv_qt_bin_dir/uic"; then
          QT_UIC="$bnv_qt_bin_dir/uic"
        else
          # Old versions of Qt don't have uic
          QT_UIC=
        fi
        QT_MOC="$bnv_qt_bin_dir/moc"
      else
      # Last possibility is that they are in $PATH
        QT_UIC="`which uic`"
        QT_MOC="`which moc`"
      fi
    fi
    # All variables are defined, report the result
    AC_MSG_RESULT([$have_qt:
    QT_CXXFLAGS=$QT_CXXFLAGS
    QT_DIR=$QT_DIR
    QT_LIBS=$QT_LIBS
    QT_UIC=$QT_UIC
    QT_MOC=$QT_MOC])
  else
    # Qt was not found
    QT_CXXFLAGS=
    QT_DIR=
    QT_LIBS=
    QT_UIC=
    QT_MOC=
    AC_MSG_RESULT($have_qt)
  fi
  AC_SUBST(QT_CXXFLAGS)
  AC_SUBST(QT_DIR)
  AC_SUBST(QT_LIBS)
  AC_SUBST(QT_UIC)
  AC_SUBST(QT_MOC)

  #### Being paranoid:
  if test x"$have_qt" = xyes; then
    AC_MSG_CHECKING(correct functioning of Qt installation)
    AC_CACHE_VAL(bnv_cv_qt_test_result,
    [
      cat > bnv_qt_test.h << EOF
#include <qobject.h>
class Test : public QObject
{
Q_OBJECT
public:
  Test() {}
  ~Test() {}
public slots:
  void receive() {}
signals:
  void send();
};
EOF

      cat > bnv_qt_main.$ac_ext << EOF
#include "bnv_qt_test.h"
#include <qapplication.h>
int main( int argc, char **argv )
{
  QApplication app( argc, argv );
  Test t;
  QObject::connect( &t, SIGNAL(send()), &t, SLOT(receive()) );
}
EOF

      bnv_cv_qt_test_result="failure"
      bnv_try_1="$QT_MOC bnv_qt_test.h -o moc_bnv_qt_test.$ac_ext >/dev/null 2>/dev/null"
      AC_TRY_EVAL(bnv_try_1)
      if test x"$ac_status" != x0; then
        echo "$bnv_err_1" >&AC_FD_CC
        echo "configure: could not run $QT_MOC on:" >&AC_FD_CC
        cat bnv_qt_test.h >&AC_FD_CC
      else
        bnv_try_2="$CXX $QT_CXXFLAGS -c $CXXFLAGS -o moc_bnv_qt_test.o moc_bnv_qt_test.$ac_ext >/dev/null 2>/dev/null"
        AC_TRY_EVAL(bnv_try_2)
        if test x"$ac_status" != x0; then
          echo "$bnv_err_2" >&AC_FD_CC
          echo "configure: could not compile:" >&AC_FD_CC
          cat moc_bnv_qt_test.$ac_ext >&AC_FD_CC
        else
          bnv_try_3="$CXX $QT_CXXFLAGS -c $CXXFLAGS -o bnv_qt_main.o bnv_qt_main.$ac_ext >/dev/null 2>/dev/null"
          AC_TRY_EVAL(bnv_try_3)
          if test x"$ac_status" != x0; then
            echo "$bnv_err_3" >&AC_FD_CC
            echo "configure: could not compile:" >&AC_FD_CC
            cat bnv_qt_main.$ac_ext >&AC_FD_CC
          else
            bnv_try_4="$CXX $QT_LIBS $LIBS -o bnv_qt_main bnv_qt_main.o moc_bnv_qt_test.o >/dev/null 2>/dev/null"
            AC_TRY_EVAL(bnv_try_4)
            if test x"$ac_status" != x0; then
              echo "$bnv_err_4" >&AC_FD_CC
            else
              bnv_cv_qt_test_result="success"
            fi
          fi
        fi
      fi
    ])dnl AC_CACHE_VAL bnv_cv_qt_test_result
    AC_MSG_RESULT([$bnv_cv_qt_test_result]);
    if test x"$bnv_cv_qt_test_result" = "xfailure"; then
      AC_MSG_ERROR([Failed to find matching components of a complete
                  Qt installation. Try using more options,
                  see ./configure --help.])
    fi

    rm -f bnv_qt_test.h moc_bnv_qt_test.$ac_ext moc_bnv_qt_test.o \
          bnv_qt_main.$ac_ext bnv_qt_main.o bnv_qt_main
  fi
])

dnl Internal subroutine of BNV_HAVE_QT
dnl Set bnv_qt_dir bnv_qt_include_dir bnv_qt_bin_dir bnv_qt_lib_dir bnv_qt_lib
AC_DEFUN([BNV_PATH_QT_DIRECT],
[
  ## Binary utilities ##
  if test x"$with_Qt_bin_dir" != x; then
    bnv_qt_bin_dir=$with_Qt_bin_dir
  fi
  ## Look for header files ##
  if test x"$with_Qt_include_dir" != x; then
    bnv_qt_include_dir="$with_Qt_include_dir"
  else
    # The following header file is expected to define QT_VERSION.
    qt_direct_test_header=qglobal.h
    # Look for the header file in a standard set of common directories.
    bnv_include_path_list="
      /usr/include
      `ls -dr ${QTDIR}/include 2>/dev/null`
      `ls -dr /usr/include/qt* 2>/dev/null`
      `ls -dr /usr/lib/qt*/include 2>/dev/null`
      `ls -dr /usr/local/qt*/include 2>/dev/null`
      `ls -dr /opt/qt*/include 2>/dev/null`
      `ls -dr /Developer/qt*/include 2>/dev/null`
    "
    for bnv_dir in $bnv_include_path_list; do
      if test -r "$bnv_dir/$qt_direct_test_header"; then
        bnv_dirs="$bnv_dirs $bnv_dir"
      fi
    done
    # Now look for the newest in this list
    bnv_prev_ver=0
    for bnv_dir in $bnv_dirs; do
      bnv_this_ver=`egrep -w '#define QT_VERSION' $bnv_dir/$qt_direct_test_header | sed s/'#define QT_VERSION'//`
      if expr $bnv_this_ver '>' $bnv_prev_ver > /dev/null; then
        bnv_qt_include_dir=$bnv_dir
        bnv_prev_ver=$bnv_this_ver
      fi
    done
  fi dnl Found header files.

  # Are these headers located in a traditional Trolltech installation?
  # That would be $bnv_qt_include_dir stripped from its last element:
  bnv_possible_qt_dir=`dirname $bnv_qt_include_dir`
  if (test -x $bnv_possible_qt_dir/bin/moc) &&
     ((ls $bnv_possible_qt_dir/lib/libqt* > /dev/null 2>/dev/null) ||
      (ls $bnv_possible_qt_dir/lib64/libqt* > /dev/null 2>/dev/null)); then
    # Then the rest is a piece of cake
    bnv_qt_dir=$bnv_possible_qt_dir
    bnv_qt_bin_dir="$bnv_qt_dir/bin"
    if test x"$with_Qt_lib_dir" != x; then
      bnv_qt_lib_dir="$with_Qt_lib_dir"
    else
      if (test -d $bnv_qt_dir/lib64); then
	bnv_qt_lib_dir="$bnv_qt_dir/lib64"
      else
	bnv_qt_lib_dir="$bnv_qt_dir/lib"
      fi
    fi
    # Only look for lib if the user did not supply it already
    if test x"$bnv_qt_lib" = xNO; then
      bnv_qt_lib="`ls $bnv_qt_lib_dir/libqt* | sed -n 1p |
                   sed s@$bnv_qt_lib_dir/lib@@ | [sed s@[.].*@@]`"
    fi
    bnv_qt_LIBS="-L$bnv_qt_lib_dir -l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
  else
    # There is no valid definition for $QTDIR as Trolltech likes to see it
    bnv_qt_dir=
    ## Look for Qt library ##
    if test x"$with_Qt_lib_dir" != x; then
      bnv_qt_lib_dir="$with_Qt_lib_dir"
      # Only look for lib if the user did not supply it already
      if test x"$bnv_qt_lib" = xNO; then
        bnv_qt_lib="`ls $bnv_qt_lib_dir/libqt* | sed -n 1p |
                     sed s@$bnv_qt_lib_dir/lib@@ | [sed s@[.].*@@]`"
      fi
      bnv_qt_LIBS="-L$bnv_qt_lib_dir -l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
    else
      # Normally, when there is no traditional Trolltech installation,
      # the library is installed in a place where the linker finds it
      # automatically.
      # If the user did not define the library name, try with qt
      if test x"$bnv_qt_lib" = xNO; then
        bnv_qt_lib=qt
      fi
      qt_direct_test_header=qapplication.h
      qt_direct_test_main="
        int argc;
        char ** argv;
        QApplication app(argc,argv);
      "
      # See if we find the library without any special options.
      # Don't add top $LIBS permanently yet
      bnv_save_LIBS="$LIBS"
      LIBS="-l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
      bnv_qt_LIBS="$LIBS"
      bnv_save_CXXFLAGS="$CXXFLAGS"
      CXXFLAGS="-I$bnv_qt_include_dir"
      AC_TRY_LINK([#include <$qt_direct_test_header>],
        $qt_direct_test_main,
      [
        # Success.
        # We can link with no special library directory.
        bnv_qt_lib_dir=
      ], [
        # That did not work. Try the multi-threaded version
        echo "Non-critical error, please neglect the above." >&AC_FD_CC
        bnv_qt_lib=qt-mt
        LIBS="-l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
        AC_TRY_LINK([#include <$qt_direct_test_header>],
          $qt_direct_test_main,
        [
          # Success.
          # We can link with no special library directory.
          bnv_qt_lib_dir=
        ], [
          # That did not work. Try the OpenGL version
          echo "Non-critical error, please neglect the above." >&AC_FD_CC
          bnv_qt_lib=qt-gl
          LIBS="-l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
          AC_TRY_LINK([#include <$qt_direct_test_header>],
            $qt_direct_test_main,
          [
            # Success.
            # We can link with no special library directory.
            bnv_qt_lib_dir=
          ], [
            # That did not work. Maybe a library version I don't know about?
            echo "Non-critical error, please neglect the above." >&AC_FD_CC
            # Look for some Qt lib in a standard set of common directories.
            bnv_dir_list="
              `echo $bnv_qt_includes | sed ss/includess`
              /lib
	      /usr/lib64
              /usr/lib
	      /usr/local/lib64
              /usr/local/lib
	      /opt/lib64
              /opt/lib
              `ls -dr /usr/lib64/qt* 2>/dev/null`
              `ls -dr /usr/lib64/qt*/lib64 2>/dev/null`
              `ls -dr /usr/lib/qt* 2>/dev/null`
              `ls -dr /usr/local/qt* 2>/dev/null`
              `ls -dr /opt/qt* 2>/dev/null`
            "
            for bnv_dir in $bnv_dir_list; do
              if ls $bnv_dir/libqt* >/dev/null 2>/dev/null; then
                # Gamble that it's the first one...
                bnv_qt_lib="`ls $bnv_dir/libqt* | sed -n 1p |
                            sed s@$bnv_dir/lib@@ | sed s/[[.]].*//`"
                bnv_qt_lib_dir="$bnv_dir"
                break
              fi
            done
            # Try with that one
            LIBS="-l$bnv_qt_lib $X_PRE_LIBS $X_LIBS -lX11 -lXext -lXmu -lXt -lXi $X_EXTRA_LIBS"
            AC_TRY_LINK([#include <$qt_direct_test_header>],
              $qt_direct_test_main,
            [
              # Success.
              # We can link with no special library directory.
              bnv_qt_lib_dir=
            ], [
              # Leave bnv_qt_lib_dir defined
            ])
          ])
        ])
      ])
      if test x"$bnv_qt_lib_dir" != x; then
        bnv_qt_LIBS="-L$bnv_qt_lib_dir $LIBS"
      else
        bnv_qt_LIBS="$LIBS"
      fi
      LIBS="$bnv_save_LIBS"
      CXXFLAGS="$bnv_save_CXXFLAGS"
    fi dnl $with_Qt_lib_dir was not given
  fi dnl Done setting up for non-traditional Trolltech installation
])

#
# EOF
#
