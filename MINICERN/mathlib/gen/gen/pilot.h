*FCA :          Fri Mar 26 17:27:50 CET 1999 by  Federico Carminati
*               define UNIX when LINUX defined


#if defined(CERNLIB_LINUX)
#ifndef CERNLIB_UNIX
#define CERNLIB_UNIX
#endif
#endif

#if (defined(CERNLIB_UNIX))&&(!defined(CERNLIB_SINGLE))
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif 

#if defined(CERNLIB_DECS)||defined(CERNLIB_QMALPH)||defined(CERNLIB_APOLLO)||defined(CERNLIB_SGI)||defined(CERNLIB_NEXT)||defined(CERNLIB_LINUX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_CONVEX32)||defined(CERNLIB_QFAPOGEE)||defined(CERNLIB_QFEPC)||defined(CERNLIB_QFMSOFT)||defined(CERNLIB_QFDEC)||defined(CERNLIB_WINNT)
#ifndef CERNLIB_NOQUAD
#define CERNLIB_NOQUAD
#endif
#endif

#if !defined(CERNLIB_NOQUAD)
#ifndef CERNLIB_QUAD
#define CERNLIB_QUAD
#endif
#endif

#if defined(CERNLIB_QMALPH)
#ifndef CERNLIB_FORTRAN
#define CERNLIB_FORTRAN
#endif
#endif

#if (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_SINGLE))
#ifndef CERNLIB_CONVEX64
#define CERNLIB_CONVEX64
#endif
#endif
#if (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_CONVEX64))
#ifndef CERNLIB_CONVEX32
#define CERNLIB_CONVEX32
#endif
#endif
#if defined(CERNLIB_CONVEX32)||defined(CERNLIB_CONVEX64)
#ifndef CERNLIB_CONVEX
#define CERNLIB_CONVEX
#endif
#endif

#if 0
* DEC Fortran 1.0 is used for Windows/NT
#endif
#if defined(CERNLIB_WINNT) && !defined(CERNLIB_QFMSOFT)
#ifndef CERNLIB_DECS
#define CERNLIB_DECS
#endif

#ifndef CERNLIB_QMALPH
#define CERNLIB_QMALPH
#endif

#ifndef CERNLIB_FORTRAN
#define CERNLIB_FORTRAN
#endif
#endif

#if (defined(CERNLIB_MSDOS))&&(!defined(CERNLIB_NDP))&&(!defined(CERNLIB_WINNT))
#ifndef CERNLIB_QF2C
#define CERNLIB_QF2C
#endif
#endif

#if defined(CERNLIB_LINUX)
#ifndef CERNLIB_QF2C
#define CERNLIB_QF2C
#endif
#endif

#if defined(CERNLIB_IBMMVS)||defined(CERNLIB_IBMVM)
#ifndef CERNLIB_IBM
#define CERNLIB_IBM
#endif
#endif

#if defined(CERNLIB_CDC)||defined(CERNLIB_CRAY)
#ifndef CERNLIB_SINGLE
#define CERNLIB_SINGLE
#endif
#endif

#if defined(CERNLIB_IBM)||defined(CERNLIB_VAX)||defined(CERNLIB_NORD)
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif

#if defined(CERNLIB_IBM)
#ifndef CERNLIB_QMIBMXA
#define CERNLIB_QMIBMXA
#endif
#endif

#if 0
* Defines for the tests not used to build lib
#endif

#if defined(CERNLIB_DOUBLE)
#ifndef CERNLIB_CMPXDOUB
#define CERNLIB_CMPXDOUB
#endif
#endif

#if defined(CERNLIB_CRAY)
#ifndef CERNLIB_CMPXDOUB
#define CERNLIB_CMPXDOUB
#endif
#endif

#if (defined(CERNLIB_UNIX))&&(!defined(CERNLIB_SINGLE))
#ifndef CERNLIB_QIEEE 
#define CERNLIB_QIEEE 
#endif
#endif

#if defined(CERNLIB_IBMRT)
#ifndef CERNLIB_IBMRS 
#define CERNLIB_IBMRS 
#endif
#endif
