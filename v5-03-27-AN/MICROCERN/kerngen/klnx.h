#if 0
* This pilot patch was created from kernlnx.car patch _klnx
* This directory was created from kernlnx.car patch qmlnx
*                 Normal Unix system machine
*                Posix call for setjmp/longjmp
*                 IEEE floating point
*       ISA standard routines, ISHFT, IOR, etc
*       Hollerith storage not orthodox
*              UCOPY et al. to copy integers
#endif
#ifndef CERNLIB_QMLNX
#define CERNLIB_QMLNX
#endif
#ifndef CERNLIB_QPOSIX
#define CERNLIB_QPOSIX
#endif
#ifndef CERNLIB_QIEEE
#define CERNLIB_QIEEE
#endif
#if (!defined(CERNLIB_PPC))
#  ifdef CERNLIB_QISASTD
#    undef CERNLIB_QISASTD
#  endif
#  ifdef CERNLIB_QORTHOLL
#    undef CERNLIB_QORTHOLL
#  endif
#else
#  ifndef CERNLIB_QISASTD
#    define CERNLIB_QISASTD
#  endif
#  ifndef CERNLIB_QORTHOLL
#    define CERNLIB_QORTHOLL
#  endif
#endif
#ifndef CERNLIB_QINTCOPY
#define CERNLIB_QINTCOPY
#endif
#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
