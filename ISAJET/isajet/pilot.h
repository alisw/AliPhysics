
#ifdef CERNLIB_CERN
#define CERNLIB_RANFCALL
#else
#define CERNLIB_RANFFTN
#define CERNLIB_NOCERN
#endif

#ifdef CERNLIB_HPUX
#ifndef CERNLIB_IMPNONE
#define CERNLIB_IMPNONE
#endif
#endif
#ifdef CERNLIB_IBMRT
#ifndef CERNLIB_IMPNONE
#define CERNLIB_IMPNONE
#endif
#endif
#ifdef CERNLIB_VAX
#ifndef CERNLIB_IMPNONE
#define CERNLIB_IMPNONE
#endif
#endif

#if !defined(CERNLIB_SINGLE)
#ifndef CERNLIB_DOUBLE
#define CERNLIB_DOUBLE
#endif
#endif

#ifndef CERNLIB_PDFLIB
#define CERNLIB_PDFLIB
#endif

#ifndef CERNLIB_STDIO
#define CERNLIB_STDIO
#endif
#ifndef CERNLIB_MOVEFTN
#define CERNLIB_MOVEFTN
#endif
