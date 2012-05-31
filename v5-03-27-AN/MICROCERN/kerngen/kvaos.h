#if 0
*      for Vax Alpha OSF, 32-bit mode
* This pilot patch was created from kernvmi.car patch _kvaos
* This directory was created from kernvmi.car patch qmvaos
* This directory was created from kernfor.car patch qmvaos
*                 Normal Unix system machine
*                 external names with underscores
*                 IEEE floating point
*               ISA standard functions available
*               MIL standard routines, IBITS, MVBITS, etc
*       Hollerith storage not orthodox
*                 Hollerith constants exist
*              EQUIVALENCE Hollerith/Character ok
*               running Unix
*               running Unix system BSD
*             use sigaction
*               BSD version for SETENVF
*               Posix call for setjmp/longjmp
#endif
#ifndef CERNLIB_QMVAO
#define CERNLIB_QMVAO
#endif
#ifndef CERNLIB_QMVAOS
#define CERNLIB_QMVAOS
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
#ifndef CERNLIB_QMILSTD
#define CERNLIB_QMILSTD
#endif
#ifdef CERNLIB_QORTHOLL
#undef CERNLIB_QORTHOLL
#endif
#ifndef CERNLIB_QS_UNIX
#define CERNLIB_QS_UNIX
#endif
#ifndef CERNLIB_QSYSBSD
#define CERNLIB_QSYSBSD
#endif
#ifndef CERNLIB_QSIGPOSIX
#define CERNLIB_QSIGPOSIX
#endif
#ifndef CERNLIB_QENVBSD
#define CERNLIB_QENVBSD
#endif
#ifndef CERNLIB_QSIGJMP
#define CERNLIB_QSIGJMP
#endif
#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
