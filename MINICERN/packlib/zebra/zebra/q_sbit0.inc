*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:24  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_SBIT0_INC
#define CERNLIB_ZEBRA_Q_SBIT0_INC
*
*
* q_sbit0.inc
*

#if defined(CERNLIB_APOLLO)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHFT(1,IZP-1)) )

#elif defined(CERNLIB_CDC)
      MSBIT0 (IZW,IZP)   = IZW .AND. .NOT.SHIFT(MASK(1),IZP)

#elif defined(CERNLIB_CRAY)
      MSBIT0 (IZW,IZP)       =  IZW .AND. SHIFT(MASK(65),IZP)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      MSBIT0 (IZW,IZP)   = IBCLR (IZW,IZP-1)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      MSBIT0 (IZW,IZP)   = IBCLR (IZW,IZP-1)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      MSBIT0 (IZZW,IZZP)   = IAND (IZZW, NOT(ISHFT(1,IZZP-1)) )

#elif defined(CERNLIB_QF_F2C)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHIFT(1,IZP-1)) )
#endif

#elif defined(CERNLIB_QF_SIEM)
      MSBIT0 (IZW,IZP)   = LAND (IZW, LCOMPL(SHFTL(1,IZP-1)) )

#elif defined(CERNLIB_LINUX)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHIFT(1,IZP-1)))

#elif defined(CERNLIB_PDP)
      MSBIT0 (IZW,IZP)   = IAND (IZW, NOT(ISHFT(1,IZP-1)) )

#elif defined(CERNLIB_SUN)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHIFT(1,IZP-1)) )

#elif defined(CERNLIB_VAXVMS)

#elif defined(CERNLIB_QMVDS)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHFT(1,IZP-1)) )

#elif 1
      MSBIT0 (IZW,IZP)   = IAND (IZW, NOT(ISHFT(1,IZP-1)) )
#endif

#endif
