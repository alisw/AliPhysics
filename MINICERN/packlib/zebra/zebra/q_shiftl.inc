*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:25  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_SHIFTL_INC
#define CERNLIB_ZEBRA_Q_SHIFTL_INC
*
*
* q_shiftl.inc
*
#if defined(CERNLIB_APOLLO)
      ISHFTL (IZW,NZB) = LSHFT (IZW, NZB)

#elif defined(CERNLIB_CDC)
      ISHFTL (IZW,NZB) = SHIFT (IZW,NZB) .AND. MASK(60-NZB)

#elif defined(CERNLIB_CRAY)
      ISHFTL (IZW,NZB) = SHIFTL (IZW, NZB)

#elif (defined(CERNLIB_MSDOS))&&(defined(CERNLIB_QF_F2C))
       ISHFTL(IZW,NZB) = LSHIFT(IZW,NZB)

#elif defined(CERNLIB_QF_SIEM)
      ISHFTL (IZW,NZB) = SHFTL (IZW, NZB)

#elif defined(CERNLIB_LINUX)
      ISHFTL (IZW,NZB) = LSHIFT (IZW,NZB)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
#elif (defined(CERNLIB_SUN))&&(defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)

#elif defined(CERNLIB_QMVDS)
C   ISHFTL / ISHFTR  are intrinsics of vdS

#elif 1
      ISHFTL (IZW,NZB) = ISHFT (IZW, NZB)

#endif

#endif
