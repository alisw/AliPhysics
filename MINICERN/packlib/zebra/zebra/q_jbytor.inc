*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:10  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_JBYTOR_INC
#define CERNLIB_ZEBRA_Q_JBYTOR_INC
*
*
* q_jbytor.inc
*
#if defined(CERNLIB_ALLIANT)
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ, IBITS(IZW,IZP-1,NZB))

#elif defined(CERNLIB_APOLLO)
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                        RSHFT (LSHFT(IZW,33-IZP-NZB), 32-NZB) )

#elif defined(CERNLIB_CDC)
      JBYTOR (MZ,IZW,IZP,NZB) = MZ .OR. (SHIFT(IZW,61-IZP)
     +                             .AND. .NOT.MASK(60-NZB) )

#elif defined(CERNLIB_CRAY)
      JBYTOR (MZ,IZW,IZP,NZB) = MZ .OR. (SHIFTR(IZW,IZP-1)
     +                             .AND. MASK(128-NZB) )

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      JBYTOR (MZ,IZW,IZP,NZB) = MZ .OR.  IBITS(IZW,IZP-1,NZB)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      JBYTOR (MZ,IZW,IZP,NZB) = MZ .OR.  IBITS(IZW,IZP-1,NZB)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      JBYTOR (MZZ,IZZW,IZZP,NZZB) = IOR  (MZZ, IBITS(IZZW,IZZP-1,NZZB))
#elif defined(CERNLIB_QF_F2C)
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                        ishftr (LSHIFT(IZW,33-IZP-NZB),32-NZB) )
#endif

#elif defined(CERNLIB_HPUX)
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ, IBITS(IZW,IZP-1,NZB))

#elif defined(CERNLIB_QF_SIEM)
      JBYTOR (MZ,IZW,IZP,NZB) = LOR  (MZ,
     +                        SHFTR (SHFTL(IZW,33-IZP-NZB), 32-NZB) )

#elif defined(CERNLIB_LINUX)
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                        ISHFTR (LSHIFT(IZW,33-IZP-NZB),32-NZB) )

#elif defined(CERNLIB_OS9)
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ, IBITS(IZW,IZP-1,NZB))

#elif defined(CERNLIB_PDP)
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ, IBITS(IZW,IZP-1,NZB))

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                       lrshft (LSHIFT(IZW,33-IZP-NZB), 32-NZB))

#elif defined(CERNLIB_QMVDS)
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ, JBYT(IZW,IZP,NZB))

#elif 1
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ,
     +                        ISHFT (ISHFT(IZW,33-IZP-NZB),-32+NZB) )
#endif

#endif
