#if 0
* This pilot patch was created from kerndos.car patch _kdos
* This directory was created from kerndos.car patch qmdos
#endif
#if (!defined(CERNLIB_QF_DEC))&&(!defined(CERNLIB_QF_NDP))
#ifndef CERNLIB_QF_F2C
#define CERNLIB_QF_F2C
#endif
#endif
#if !defined(CERNLIB_QS_WNT)
#ifndef CERNLIB_QS_DOS
#define CERNLIB_QS_DOS
#endif
#endif
#ifndef CERNLIB_QMDOS
#define CERNLIB_QMDOS
#endif
#ifndef CERNLIB_QX_SC
#define CERNLIB_QX_SC
#endif
#ifndef CERNLIB_QIEEE
#define CERNLIB_QIEEE
#endif
#if defined(CERNLIB_QF_NDP)||defined(CERNLIB_QF_DEC)
#ifndef CERNLIB_QISASTD
#define CERNLIB_QISASTD
#endif
#ifndef CERNLIB_QMILSTD
#define CERNLIB_QMILSTD
#endif
#endif
#ifdef CERNLIB_QORTHOLL
#undef CERNLIB_QORTHOLL
#endif
#ifndef CERNLIB_F77TRARG
#define CERNLIB_F77TRARG
#endif
#ifndef CERNLIB_QCFIO
#define CERNLIB_QCFIO
#endif
#ifndef CERNLIB_QGETCWD
#define CERNLIB_QGETCWD
#endif
#ifndef CERNLIB_QINTCOPY
#define CERNLIB_QINTCOPY
#endif
#ifndef CERNLIB_QINTZERO
#define CERNLIB_QINTZERO
#endif
