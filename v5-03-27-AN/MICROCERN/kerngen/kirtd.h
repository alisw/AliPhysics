#if 0
* This directory was created from kernirt.car patch qmirtd
* This pilot patch was created from kernirt.car patch _kirtd
*                 Normal Unix system machine
*               Posix call for setjmp/longjmp
*               MIL standard routines, IBITS, MVBITS, etc
*               ISA standard routines, ISHFT, IOR, etc
*                 External names with underscore
*                 IEEE floating point
*                 Hollerith constants exist
*              EQUIVALENCE Hollerith/Character ok
*              Orthodox Hollerith storage left to right
#endif
#ifndef CERNLIB_QMIRTD
#define CERNLIB_QMIRTD
#endif
#ifndef CERNLIB_QSIGJMP
#define CERNLIB_QSIGJMP
#endif
#ifndef CERNLIB_QMILSTD
#define CERNLIB_QMILSTD
#endif
#ifndef CERNLIB_QISASTD
#define CERNLIB_QISASTD
#endif
#ifndef CERNLIB_QX_SC
#define CERNLIB_QX_SC
#endif
#ifndef CERNLIB_QIEEE
#define CERNLIB_QIEEE
#endif
#ifndef CERNLIB_QORTHOLL
#define CERNLIB_QORTHOLL
#endif
#ifdef CERNLIB_QINTCOPY
#undef CERNLIB_QINTCOPY
#endif
#ifdef CERNLIB_QINTZERO
#undef CERNLIB_QINTZERO
#endif
