*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:15:06  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_JBIT_INC
#define CERNLIB_ZEBRA_Q_JBIT_INC
*
*
* q_jbit.inc
*
#if defined(CERNLIB_ALLIANT)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_APOLLO)
      JBIT (IZW,IZP)    =   AND (RSHFT(IZW,IZP-1), 1)

#elif defined(CERNLIB_CDC)
      JBIT (IZW,IZP)     = SHIFT(IZW,61-IZP) .AND. .NOT.MASK(59)

#elif defined(CERNLIB_CRAY)
      JBIT (IZW,IZP)     = SHIFTR(IZW,IZP-1) .AND. MASK(127)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      JBIT (IZZW,IZZP)      = IBITS (IZZW,IZZP-1,1)
#elif defined(CERNLIB_QF_F2C)
      JBIT (IZW,IZP)      = AND(RSHIFT(IZW,IZP-1),1)
#endif

#elif defined(CERNLIB_HPUX)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_QF_SIEM)
      JBIT (IZW,IZP)    = SHFTR (SHFTL(IZW,32-IZP), 31)

#elif defined(CERNLIB_LINUX)
      JBIT(IZW,IZP)     = AND(ISHFTR(IZW,IZP-1),1)

#elif defined(CERNLIB_OS9)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_PDP)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      JBIT (IZW,IZP)    = AND (lrshft(IZW,IZP-1), 1)

#elif defined(CERNLIB_QMUNO)
      JBIT (IZW,IZP)     = BITS (IZW,37-IZP,1)

#elif defined(CERNLIB_QMVDS)
      JBIT (IZW,IZP)  =  JBYT (IZW,IZP,1)

#elif defined(CERNLIB_VAXVMS)

#elif 1
      JBIT (IZW,IZP)     =  IAND (ISHFT(IZW,1-IZP), 1)

#endif

#endif
