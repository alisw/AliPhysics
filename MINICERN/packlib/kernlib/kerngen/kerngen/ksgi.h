#if 0
*       for Silicon Graphics Iris
* This pilot patch was created from kernsgi.car patch _ksgi
* This directory was created from kernsgi.car patch qmsgi
* This directory was created from kernfor.car patch qmsgi
*        there are still problems with IBITS
*                 Normal Unix system machine
*                 external names with underscore
*                 floating point is IEEE
*               ISA standard routines, ISHFT, IOR, etc
*               MIL standard routines, IBITS, MVBITS, etc
*                 Hollerith constants exist
*              EQUIVALENCE Hollerith/Character ok
*              Orthodox Hollerith storage left to right
*              Internal double-precision
*               running Unix
*               Posix call for setjmp/longjmp
#endif
#ifdef CERNLIB_QMILSTD
#undef CERNLIB_QMILSTD
#endif
#ifndef CERNLIB_QMSGI
#define CERNLIB_QMSGI
#endif
#ifndef CERNLIB_QX_SC
#define CERNLIB_QX_SC
#endif
#ifndef CERNLIB_QIEEE
#define CERNLIB_QIEEE
#endif
#ifndef CERNLIB_QISASTD
#define CERNLIB_QISASTD
#endif
#ifndef CERNLIB_QORTHOLL
#define CERNLIB_QORTHOLL
#endif
#ifndef CERNLIB_INTDOUBL
#define CERNLIB_INTDOUBL
#endif
#ifndef CERNLIB_QS_UNIX
#define CERNLIB_QS_UNIX
#endif
#ifndef CERNLIB_QSIGJMP
#define CERNLIB_QSIGJMP
#endif
#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
