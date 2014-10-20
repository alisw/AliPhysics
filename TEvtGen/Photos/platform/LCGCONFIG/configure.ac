#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT( [Photos C++ Interface], [3.3], [tomasz.przedzinski@cern.ch] )

AC_CONFIG_AUX_DIR([config])

AM_INIT_AUTOMAKE

AC_LANG_CPLUSPLUS

AC_PREFIX_DEFAULT([`pwd`])
AC_CONFIG_HEADERS([config/config.h])
AC_CONFIG_SRCDIR([src/photosCInterfaces/Photos.cxx])
##############################################################################
# Check for debug flags

AC_MSG_CHECKING([if debug flags are enabled])
AC_ARG_ENABLE(
  [debug],
  [AS_HELP_STRING([--enable-debug],[Enable debug flags])],
  [
    AC_MSG_RESULT([$enable_debug]) ;
    FFLAGS="-g3 -O0" ;
    CXXFLAGS="-g3 -O0"
  ],
  [
    AC_MSG_RESULT([no]) ;
    FFLAGS="-O2" ;
    CXXFLAGS="-O2"
  ]
)

# Check for HepMC

AC_MSG_CHECKING([if HepMC location defined])
AC_ARG_WITH(
  [hepmc],
  [AS_HELP_STRING([--with-hepmc], [Path to HepMC 2.xx (use --without-hepmc to compile without HepMC)] )],
  [
    AC_MSG_RESULT([$with_hepmc]) ;
    AC_SUBST(HEPMC_DIR,[$with_hepmc])
  ],
  [AC_MSG_RESULT([no])]
)

# autoconf 2.59 requires at least one check to be preformed unconditionally
#               before other conditional checks can be done
AC_CHECK_HEADERS([stdio.h],,)

if ! test "$with_hepmc" == "no"; then

CPPFLAGS="-I$with_hepmc/include $CPPFLAGS"

AC_CHECK_HEADERS([HepMC/Version.h], ,
[ echo "" &&
  echo "ERROR: HepMC not found. To provide location use:" &&
  echo "       './configure --with-hepmc=<path_to_HepMC>' "
  echo "       or use:"
  echo "       './configure --without-hepmc"
  echo "       (Consider using 'source platform/afs.paths.sh' (or .csh) )"
  echo ""
  exit -1
])

else
  with_hepmc=""
fi

##############################################################################
# Check for programs
AC_PROG_MAKE_SET
AC_PROG_LN_S
AC_PROG_CC
AC_PROG_CXX

if ${CXX} --version | grep -E "GCC.* 3" >&/dev/null ; then
	AC_PROG_F77([g77 gfortran f77 fortran fort])
else
	AC_PROG_F77([gfortran g77 f77 fortran fort])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getcwd pow sqrt], ,
[ echo "" &&
  echo "ERROR: Basic functions not found. Something's wrong with paths." &&
  exit -1
])

AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_LN_S

##############################################################################
# Check for other paths

AC_MSG_CHECKING([if Pythia8 location defined])
AC_ARG_WITH(
  [pythia8],
  [AS_HELP_STRING([--with-pythia8], [Path to Pythia 8.1 or later] )],
  [
    AC_MSG_RESULT([$with_pythia8]) ;
    AC_SUBST(PYTHIA8_DIR,[$with_pythia8])
  ],
  [
    AC_MSG_RESULT([no]) ;
    echo ""
    echo "INFO: Pythia8 location not specified. Not all examples will be available." ;
    echo "      To define location use parameter:  --with-pythia8" ;
    echo ""
  ]
)

if test "x$with_pythia8" != "x"; then
  if test "x$PYTHIA8DATA" == "x" ; then
    echo ""
    echo "WARNING: Pythia8 location specified but PYTHIA8DATA variable not set."
    echo "         This variable is required to run the examples. Try:"
    echo "         export PYTHIA8DATA=<path_to_pythia8>/xmldoc"
    echo ""
  fi
fi

AC_MSG_CHECKING([if MC-TESTER location defined])
AC_ARG_WITH(
  [mc-tester],
  [AS_HELP_STRING([--with-mc-tester], [Path to MC-TESTER 1.23 or later] )],
  [
    AC_MSG_RESULT([$with_mc_tester]) ;
    AC_SUBST(MCTESTER_DIR,[$with_mc_tester])
  ],
  [
    AC_MSG_RESULT([no]) ;
    echo "" ;
    echo "INFO: MC-TESTER location not specified. Not all examples will be available." ;
    echo "      To define location use parameter:  --with-mc-tester" ;
    echo ""
  ]
)

if test "x$with_mc_tester" != "x"; then
  AC_CHECK_PROG(HAS_ROOT_CONFIG,root-config,[yes],[no])
  if test "$HAS_ROOT_CONFIG" == "no"; then
    echo ""
    echo "WARNING: root-config not found. ROOT is not installed or path"
    echo "         to ROOT bin folder not provied via \$PATH variable. Use:"
    echo "         export PATH=\$PATH:<path_to_root>/bin"
    echo "         export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\`root-config --libdir\`"
    echo ""
  else
    ROOT_LIBS=`root-config --glibs`
    ROOT_INC=`root-config --incdir`
    AC_SUBST(ROOTLIBS,[$ROOT_LIBS])
    AC_SUBST(ROOTINC,[$ROOT_INC])
  fi
fi

AC_MSG_CHECKING([if Tauola C++ Interface location defined])
AC_ARG_WITH(
  [tauola],
  [AS_HELP_STRING([--with-tauola], [Path to Tauola C++ Interface v1.0.0 or later] )],
  [
    AC_MSG_RESULT([$with_tauola]) ;
    AC_SUBST(TAUOLA_DIR,[$with_tauola])
  ],
  [
    AC_MSG_RESULT([no]) ;
    echo ""
    echo "INFO: Tauola location not specified. Not all examples will be available." ;
    echo "      To define location use parameter:  --with-tauola" ;
    echo ""
  ]
)

AM_CONDITIONAL([HAS_HEPMC],    [test "x$with_hepmc"     != "x"])
AM_CONDITIONAL([HAS_PYTHIA8],  [test "x$with_pythia8"   != "x"])
AM_CONDITIONAL([HAS_MCTESTER], [test "x$with_mc_tester" != "x"])
AM_CONDITIONAL([HAS_TAUOLA],   [test "x$with_tauola"    != "x"])

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])
AC_OUTPUT
