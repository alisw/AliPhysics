#if 0
*-- Author :
#endif
#if defined(CERNLIB_MSDOS)||defined(CERNLIB_LINUX)
#ifndef CERNLIB_F2C
#define CERNLIB_F2C
#endif
#endif
#if defined(CERNLIB_SUN)||defined(CERNLIB_SGI)||defined(CERNLIB_DECS)||defined(CERNLIB_CONVEX)||defined(CERNLIB_IBMRT)||defined(CERNLIB_AIX370)
#ifndef CERNLIB_UNIX
#define CERNLIB_UNIX
#endif
#endif
#if defined(CERNLIB_HPUX)||defined(CERNLIB_APOLLO)||defined(CERNLIB_IPSC)||defined(CERNLIB_NEXT)
#ifndef CERNLIB_UNIX
#define CERNLIB_UNIX
#endif
#endif
#if defined(CERNLIB_IBM)||defined(CERNLIB_IBMMVS)||defined(CERNLIB_AIX370)
#ifndef CERNLIB_IBMALL
#define CERNLIB_IBMALL
#endif
#endif
#if defined(CERNLIB_APOLLO)||defined(CERNLIB_IBMALL)||defined(CERNLIB_VAX)
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif

#ifdef CERNLIB_WINNT
#  ifndef CERNLIB_UNIX
#    define CERNLIB_UNIX
#  endif
#  ifndef CERNLIB_DOUBLE
#    define CERNLIB_DOUBLE
#  endif
#endif

#if (defined(CERNLIB_UNIX))&&(!defined(CERNLIB_SINGLE))
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif

#if defined(CERNLIB_UNIX)||defined(CERNLIB_QMALPH)
#if !defined(CERNLIB_BSLASH) && !defined(CERNLIB_QFMSOFT)
#define CERNLIB_BSLASH
#endif
#endif

#if defined(CERNLIB_UNIX)
#ifndef CERNLIB_USRJMP
#define CERNLIB_USRJMP
#endif
#endif

#if defined(CERNLIB_BLDLIB)
#ifndef CERNLIB_HIGZ
#define CERNLIB_HIGZ
#endif
#ifndef CERNLIB_CG
#define CERNLIB_CG
#endif

#ifndef CERNLIB_MONITOR
#define CERNLIB_MONITOR
#endif

#ifndef CERNLIB_FLUKA
#define CERNLIB_FLUKA
#endif
#ifndef CERNLIB_COMIS
#define CERNLIB_COMIS
#endif
#ifndef CERNLIB_DZDOC
#define CERNLIB_DZDOC
#endif
#endif
