*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:25  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_SHIFTR_INC
#define CERNLIB_ZEBRA_Q_SHIFTR_INC
*
*
* q_shiftr.inc
*
#if defined(CERNLIB_APOLLO)
      ISHFTR (IZW,NZB) = RSHFT (IZW, NZB)

#elif defined(CERNLIB_CDC)
      ISHFTR (IZW,NZB) = SHIFT (IZW.AND.MASK(60-NZB), 60-NZB)

#elif defined(CERNLIB_CRAY)
      ISHFTR (IZW,NZB) = SHIFTR (IZW, NZB)

#elif (defined(CERNLIB_MSDOS))&&(defined(CERNLIB_QF_F2C))
C--
C--    Note: RSHIFT is arithmetic with F2C, it cannot be used
C
C       ISHFTR(IZW,NZB) = RSHIFT(IZW,NZB)
C       ISHFT(IZW,NZB)  = LSHIFT(IZW,NZB)

#elif defined(CERNLIB_QF_SIEM)
      ISHFTR (IZW,NZB) = SHFTR (IZW, NZB)

#elif defined(CERNLIB_LINUX)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      ISHFTR (IZW,NZB) = lrshft (IZW, NZB)
#elif (defined(CERNLIB_SUN))&&(defined(CERNLIB_BUGLRSHFT))
      ISHFTR (IZW,NZB) =  ishft (IZW,-NZB)

#elif defined(CERNLIB_QMVDS)
C   ISHFTL / ISHFTR  are intrinsics of vdS

#elif 1
      ISHFTR (IZW,NZB) = ISHFT (IZW,-NZB)

#endif

#endif
