*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:15:07  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_JBYTET_INC
#define CERNLIB_ZEBRA_Q_JBYTET_INC
*
*
* q_jbytet.inc
*
#if defined(CERNLIB_ALLIANT)
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ, IBITS(IZW,IZP-1,NZB))

#elif defined(CERNLIB_APOLLO)
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                        RSHFT (LSHFT(IZW,33-IZP-NZB), 32-NZB) )

#elif defined(CERNLIB_CDC)
      JBYTET (MZ,IZW,IZP,NZB) = MZ .AND. SHIFT(IZW,61-IZP)
     +                             .AND. .NOT.MASK(60-NZB)

#elif defined(CERNLIB_CRAY)
      JBYTET (MZ,IZW,IZP,NZB) = MZ .AND. SHIFTR(IZW,IZP-1)
     +                             .AND. MASK(128-NZB)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      JBYTET (MZ,IZW,IZP,NZB) = MZ .AND. IBITS(IZW,IZP-1,NZB)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      JBYTET (MZ,IZW,IZP,NZB) = MZ .AND. IBITS(IZW,IZP-1,NZB)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      JBYTET (MZZ,IZZW,IZZP,NZZB) = IAND (MZZ, IBITS(IZZW,IZZP-1,NZZB))
#elif defined(CERNLIB_QF_F2C)
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                        ishftr (LSHIFT(IZW,33-IZP-NZB),32-NZB) )
#endif

#elif defined(CERNLIB_HPUX)
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ, IBITS(IZW,IZP-1,NZB))

#elif defined(CERNLIB_QF_SIEM)
      JBYTET (MZ,IZW,IZP,NZB) = LAND (MZ,
     +                        SHFTR (SHFTL(IZW,33-IZP-NZB), 32-NZB) )

#elif defined(CERNLIB_LINUX)
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                        ISHFTR (LSHIFT(IZW,33-IZP-NZB),32-NZB) )

#elif defined(CERNLIB_OS9)
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ, IBITS(IZW,IZP-1,NZB))

#elif defined(CERNLIB_PDP)
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ, IBITS(IZW,IZP-1,NZB))

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                       lrshft (LSHIFT(IZW,33-IZP-NZB), 32-NZB))

#elif defined(CERNLIB_QMVDS)
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ, JBYT(IZW,IZP,NZB))

#elif 1
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ,
     +                        ISHFT (ISHFT(IZW,33-IZP-NZB),-32+NZB) )
#endif

#endif
