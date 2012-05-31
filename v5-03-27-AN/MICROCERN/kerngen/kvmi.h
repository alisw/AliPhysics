#if 0
*       for Vax / Mips
* This pilot patch was created from kernvmi.car patch _kvmi
* This directory was created from kernvmi.car patch qmvmi
* This directory was created from kernfor.car patch qmvmi
*                 Normal Unix system machine
*                 external names with underscores
*                 IEEE floating point
*               ISA standard functions available
*       Hollerith storage not orthodox
*                 Hollerith constants exist
*              EQUIVALENCE Hollerith/Character ok
*               running Unix
*               running Unix system BSD
*               signal handling with BSD sigvec
*               BSD version for SETENVF
*               Posix call for setjmp/longjmp
#endif
#ifndef CERNLIB_QMVMI
#define CERNLIB_QMVMI
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
#ifdef CERNLIB_QORTHOLL
#undef CERNLIB_QORTHOLL
#endif
#ifndef CERNLIB_QS_UNIX
#define CERNLIB_QS_UNIX
#endif
#ifndef CERNLIB_QSYSBSD
#define CERNLIB_QSYSBSD
#endif
#ifndef CERNLIB_QSIGBSD
#define CERNLIB_QSIGBSD
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
