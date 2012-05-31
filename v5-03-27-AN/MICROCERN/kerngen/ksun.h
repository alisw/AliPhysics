#if 0
*       for SUN
* This pilot patch was created from kernsun.car patch _ksun
* This directory was created from kernsun.car patch qmsun
* This directory was created from kernfor.car patch qmsun
*                 Normal Unix system machine
*                 IEEE floating point
*                 external names with underscore
*                 Hollerith constants exist
*              EQUIVALENCE Hollerith/Character ok
*              Orthodox Hollerith storage left to right
*              Internal double-precision
*             signal handling with Posix sigaction
*               running Unix
#endif

#if !defined(CERNLIB_SOLARIS)
#if 0
CERNLIB_BUGLRSHFT to get round the lrshft bug in Sun f77 3.0.x
#endif
#define CERNLIB_BUGLRSHFT
#endif

#ifndef CERNLIB_QMSUN
#define CERNLIB_QMSUN
#endif

#ifndef CERNLIB_QIEEE
#define CERNLIB_QIEEE
#endif

#ifndef CERNLIB_QX_SC
#define CERNLIB_QX_SC
#endif

#ifndef CERNLIB_QORTHOLL
#define CERNLIB_QORTHOLL
#endif

#ifndef CERNLIB_INTDOUBL
#define CERNLIB_INTDOUBL
#endif

#ifndef CERNLIB_QSIGPOSIX
#define CERNLIB_QSIGPOSIX
#endif

#if defined(CERNLIB_SOLARIS)

#ifndef CERNLIB_QSIGJMP
#define CERNLIB_QSIGJMP
#endif

#ifndef CERNLIB_QGETCWD
#define CERNLIB_QGETCWD
#endif

#ifndef CERNLIB_QSIGPOSIX
#define CERNLIB_QSIGPOSIX
#endif

#ifdef CERNLIB_QSYSBSD
#undef CERNLIB_QSYSBSD
#endif

#ifdef CERNLIB_QENVBSD
#undef CERNLIB_QENVBSD
#endif

#endif

#ifndef CERNLIB_QS_UNIX
#define CERNLIB_QS_UNIX
#endif

#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
