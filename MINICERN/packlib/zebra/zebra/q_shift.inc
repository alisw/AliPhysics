*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_SHIFT_INC
#define CERNLIB_ZEBRA_Q_SHIFT_INC
*
*
* q_shift.inc
*
#if defined(CERNLIB_APOLLO)
      ISHFTL (IZW,NZB) = LSHFT (IZW, NZB)
      ISHFTR (IZW,NZB) = RSHFT (IZW, NZB)

#elif defined(CERNLIB_CDC)
      ISHFTL (IZW,NZB) = SHIFT (IZW,NZB) .AND. MASK(60-NZB)
      ISHFTR (IZW,NZB) = SHIFT (IZW.AND.MASK(60-NZB), 60-NZB)

#elif defined(CERNLIB_CRAY)
      ISHFTL (IZW,NZB) = SHIFTL (IZW, NZB)
      ISHFTR (IZW,NZB) = SHIFTR (IZW, NZB)

#elif (defined(CERNLIB_MSDOS))&&(defined(CERNLIB_QF_F2C))
C--
C--    Note: RSHIFT is arithmetic with F2C, it cannot be used
C
       ISHFTL(IZW,NZB) = LSHIFT(IZW,NZB)
C       ISHFTR(IZW,NZB) = RSHIFT(IZW,NZB)
C       ISHFT(IZW,NZB)  = LSHIFT(IZW,NZB)

#elif defined(CERNLIB_QF_SIEM)
      ISHFTL (IZW,NZB) = SHFTL (IZW, NZB)
      ISHFTR (IZW,NZB) = SHFTR (IZW, NZB)

#elif defined(CERNLIB_LINUX)
      ISHFTL (IZW,NZB) = LSHIFT (IZW,NZB)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
      ISHFTR (IZW,NZB) = lrshft (IZW, NZB)
#elif (defined(CERNLIB_SUN))&&(defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
      ISHFTR (IZW,NZB) =  ishft (IZW,-NZB)

#elif defined(CERNLIB_QMVDS)
C   ISHFTL / ISHFTR  are intrinsics of vdS

#elif 1
      ISHFTL (IZW,NZB) = ISHFT (IZW, NZB)
      ISHFTR (IZW,NZB) = ISHFT (IZW,-NZB)

#endif

#endif
