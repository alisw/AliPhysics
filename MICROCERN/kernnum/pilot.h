#if defined(CERNLIB_IBMRT)
#ifndef CERNLIB_IBMRS
#define CERNLIB_IBMRS
#endif
#endif
#if defined(CERNLIB_VAXVMS)||defined(CERNLIB_VAXULTRIX)
#ifndef CERNLIB_VAX
#define CERNLIB_VAX
#endif
#endif
#if defined(CERNLIB_VAX)
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif
#if (defined(CERNLIB_UNIX))&&(!defined(CERNLIB_SINGLE))
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif

#if defined(CERNLIB_LINUX)||defined(CERNLIB_MSDOS) && !defined(CERNLIB_WINNT)
#ifndef CERNLIB_NUMIB2
#define CERNLIB_NUMIB2
#endif
#ifdef CERNLIB_NUMD38
#undef CERNLIB_NUMD38
#endif
#ifndef CERNLIB_NUMD279
#define CERNLIB_NUMD279
#endif
#endif

#if defined(CERNLIB_VAXVMS)
#ifndef CERNLIB_NUMDE
#define CERNLIB_NUMDE
#endif
#ifndef CERNLIB_NUMD38
#define CERNLIB_NUMD38
#endif
#ifdef CERNLIB_NUMD279
#undef CERNLIB_NUMD279
#endif
#endif
#if (defined(CERNLIB_UNIX)) || defined(CERNLIB_WINNT) && (!defined(CERNLIB_QF2C))
#ifndef CERNLIB_NUMAP
#define CERNLIB_NUMAP
#endif
#ifndef CERNLIB_NUMD38
#define CERNLIB_NUMD38
#endif
#ifdef CERNLIB_NUMD279
#undef CERNLIB_NUMD279
#endif
#endif
#ifndef CERNLIB_NUMLOPRE
#define CERNLIB_NUMLOPRE
#endif
#ifdef CERNLIB_NUMHIPRE
#undef CERNLIB_NUMHIPRE
#endif
#ifndef CERNLIB_NUMRDBLE
#define CERNLIB_NUMRDBLE
#endif
#ifndef CERNLIB_NUMCDBLE
#define CERNLIB_NUMCDBLE
#endif
#ifndef CERNLIB_NUME38
#define CERNLIB_NUME38
#endif
#ifdef CERNLIB_NUME75
#undef CERNLIB_NUME75
#endif
#ifdef CERNLIB_NUME293
#undef CERNLIB_NUME293
#endif
#ifdef CERNLIB_NUME2465
#undef CERNLIB_NUME2465
#endif
#ifdef CERNLIB_NUMD75
#undef CERNLIB_NUMD75
#endif
#ifdef CERNLIB_NUMD2465
#undef CERNLIB_NUMD2465
#endif
