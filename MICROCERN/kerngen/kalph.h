#if 0
*       for VAX
* This pilot patch was created from kernvax.car patch _kalph
* This directory was created from kernvax.car patch vaxalpha
* This directory was created from kernvax.car patch qmvax
* This directory was created from kernfor.car patch qmvax
*               CC assumed available
*     external names without underscores
*     ISA standard routines, ISHFT, IOR, etc
*     MIL standard routines, IBITS, MVBITS, ISHFTC
*       Hollerith constants exist
*    EQUIVALENCE Hollerith/Character ok
*       Hollerith storage not orthodox
#endif
#ifndef CERNLIB_QMALPH
#define CERNLIB_QMALPH
#endif
#ifndef CERNLIB_QMVAXCC
#define CERNLIB_QMVAXCC
#endif
#ifndef CERNLIB_QMVAX
#define CERNLIB_QMVAX
#endif
#ifndef CERNLIB_QXNO_SC
#define CERNLIB_QXNO_SC
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
#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
