*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:07  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_JBYT_INC
#define CERNLIB_ZEBRA_Q_JBYT_INC
*
*
* q_jbyt.inc
*
#if defined(CERNLIB_ALLIANT)
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)

#elif defined(CERNLIB_APOLLO)
      JBYT (IZW,IZP,NZB)= RSHFT (LSHFT(IZW,33-IZP-NZB), 32-NZB)

#elif defined(CERNLIB_CDC)
      JBYT (IZW,IZP,NZB) = SHIFT(IZW,61-IZP) .AND. .NOT.MASK(60-NZB)

#elif defined(CERNLIB_CRAY)
      JBYT (IZW,IZP,NZB) = SHIFTR(IZW,IZP-1) .AND. MASK(128-NZB)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      JBYT (IZZW,IZZP,NZZB) = IBITS (IZZW,IZZP-1,NZZB)
#elif defined(CERNLIB_QF_F2C)
      JBYT (IZW,IZP,NZB) =  ISHFTR(LSHIFT(IZW,33-IZP-NZB),32-NZB)
#endif

#elif defined(CERNLIB_HPUX)
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)

#elif defined(CERNLIB_QF_SIEM)
      JBYT (IZW,IZP,NZB)= SHFTR (SHFTL(IZW,33-IZP-NZB), 32-NZB)

#elif defined(CERNLIB_LINUX)
      JBYT(IZW,IZP,NZB) = ISHFTR(LSHIFT(IZW,33-IZP-NZB),32-NZB)

#elif defined(CERNLIB_OS9)
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)

#elif defined(CERNLIB_PDP)
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      JBYT (IZW,IZP,NZB)= lrshft (LSHIFT(IZW,33-IZP-NZB),32-NZB)

#elif defined(CERNLIB_QMUNO)
      JBYT (IZW,IZP,NZB) = BITS (IZW,38-IZP-NZB,NZB)

#elif defined(CERNLIB_QMVDS)

#elif defined(CERNLIB_VAXVMS)

#elif 1
      JBYT (IZW,IZP,NZB) = ISHFT (ISHFT(IZW,33-IZP-NZB), -32+NZB)

#endif

#endif
