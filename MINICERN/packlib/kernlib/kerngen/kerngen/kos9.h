#if 0
*               for Microware OS-9
* This pilot patch was created from kernos9.car patch _kos9
* This directory was created from kernos9.car patch qmos9
* This directory was created from kernfor.car patch qmos9
*    external names with underscore
*                 IEEE floating point
*                 Hollerith constants exist
*              EQUIVALENCE Hollerith/Character ok
*              Internal double-precision
*              Orthodox Hollerith storage left to right
*               ISA standard routines, ISHFT, IOR, etc
*               MIL standard routines, IBITS, MVBITS, etc
*               running Unix
#endif
#ifndef CERNLIB_QMOS9
#define CERNLIB_QMOS9
#endif
#if !defined(CERNLIB_QXNO_SC)
#ifndef CERNLIB_QX_SC
#define CERNLIB_QX_SC
#endif
#endif
#ifndef CERNLIB_QIEEE
#define CERNLIB_QIEEE
#endif
#ifndef CERNLIB_QORTHOLL
#define CERNLIB_QORTHOLL
#endif
#ifndef CERNLIB_INTDOUBL
#define CERNLIB_INTDOUBL
#endif
#ifndef CERNLIB_QISASTD
#define CERNLIB_QISASTD
#endif
#ifndef CERNLIB_QMILSTD
#define CERNLIB_QMILSTD
#endif
#ifdef CERNLIB_NOGETWD
#undef CERNLIB_NOGETWD
#endif
#ifndef CERNLIB_QS_UNIX
#define CERNLIB_QS_UNIX
#endif
#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
