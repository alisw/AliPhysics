#ifdef WIN32
# ifdef CERNLIB_MSSTDCALL
#  define F77_UCASE
#  define type_of_call _stdcall
#  ifndef CERNLIB_QXCAPT
#    define CERNLIB_QXCAPT
#  endif
# else
#  define F77_LCASE
#  ifndef CERNLIB_QXNO_SC
#    define CERNLIB_QXNO_SC
#  endif
# endif
#endif

#ifndef type_of_call
#define type_of_call
#endif

