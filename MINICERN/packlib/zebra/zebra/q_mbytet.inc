*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:11  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_MBYTET_INC
#define CERNLIB_ZEBRA_Q_MBYTET_INC
*
*
* q_mbytet.inc
*
#if defined(CERNLIB_ALLIANT)
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW,
     +           NOT(ishft (ISHFT(NOT(MZ),32-NZB),-33+IZP+NZB)) )

#elif defined(CERNLIB_APOLLO)
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW, NOT(RSHFT (
     +                          LSHFT(NOT(MZ),32-NZB), 33-IZP-NZB)) )

#elif defined(CERNLIB_CDC)
      MBYTET (MZ,IZW,IZP,NZB) = IZW .AND.
     +                          SHIFT (MZ .OR. MASK(60-NZB),IZP-1)

#elif defined(CERNLIB_CRAY)
      MBYTET (MZ,IZW,IZP,NZB) = IZW .AND.
     +                             SHIFT (MZ.OR. MASK(64-NZB), IZP-1)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      MBYTET (MZZ,IZZW,IZZP,NZZB) = IAND (IZZW,
     +           NOT(ishft (ISHFT(NOT(MZZ),32-NZZB),-33+IZZP+NZZB)) )
#elif defined(CERNLIB_QF_F2C)
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW,
     +           NOT(ishftr (LSHIFT(NOT(MZ),32-NZB),33-IZP-NZB)) )
#endif

#elif defined(CERNLIB_HPUX)
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW,
     +           NOT(ishft (ISHFT(NOT(MZ),32-NZB),-33+IZP+NZB)) )

#elif defined(CERNLIB_QF_SIEM)
      MBYTET (MZ,IZW,IZP,NZB) = LAND (IZW, LCOMPL(SHFTR (
     +                       SHFTL(LCOMPL(MZ),32-NZB), 33-IZP-NZB)) )

#elif defined(CERNLIB_LINUX)
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW,
     +           NOT(ISHFTR (LSHIFT(NOT(MZ),32-NZB),33-IZP-NZB)) )

#elif defined(CERNLIB_OS9)
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW,
     +           NOT(ishft (ISHFT(NOT(MZ),32-NZB),-33+IZP+NZB)) )

#elif defined(CERNLIB_PDP)
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW, NOT(ISHFT (
     +                          ISHFT(NOT(MZ),36-NZB),-37+IZP+NZB)) )

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW, NOT(lrshft (
     +                          LSHIFT(NOT(MZ),32-NZB), 33-IZP-NZB)))

#elif defined(CERNLIB_QMVDS)
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW, NOT(RSHFT (
     +                          LSHFT(NOT(MZ),32-NZB), 33-IZP-NZB)) )

#elif 1
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW,
     +           NOT(ISHFT (ISHFT(NOT(MZ),32-NZB),-33+IZP+NZB)) )
#endif

#endif
