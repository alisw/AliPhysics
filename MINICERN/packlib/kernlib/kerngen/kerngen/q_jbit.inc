*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:35:39  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:29  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:19  mclareni
* Kernlib
*
*
*
* q_jbit.inc
*
#if defined(CERNLIB_QMLNX)
      JBIT(IZW,IZP)     = AND(ISHFTR(IZW,IZP-1),1)
      JBYT(IZW,IZP,NZB) = ISHFTR(LSHIFT(IZW,33-IZP-NZB),32-NZB)
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C)
      JBIT (IZW,IZP)      = AND(RSHIFT(IZW,IZP-1),1)
      JBYT (IZW,IZP,NZB) =  ISHFTR(LSHIFT(IZW,33-IZP-NZB),32-NZB)
#elif defined(CERNLIB_QMSUN)&&defined(CERNLIB_BUGLRSHFT)
#elif defined(CERNLIB_QMSUN)&&(!defined(CERNLIB_BUGLRSHFT))
      JBIT (IZW,IZP)    = AND (lrshft(IZW,IZP-1), 1)
      JBYT (IZW,IZP,NZB)= lrshft (LSHIFT(IZW,33-IZP-NZB),32-NZB)
#elif defined(CERNLIB_QMOS9)||defined(CERNLIB_QMHPX)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)
#elif defined(CERNLIB_QMDOS)&&(!defined(CERNLIB_QF_F2C)) || defined(CERNLIB_QFMSOFT)
      JBIT (IZW,IZP)     =  IAND (ISHFT(IZW,1-IZP), 1)
      JBYT (IZW,IZP,NZB) = ISHFT (ISHFT(IZW,33-IZP-NZB), -32+NZB)
#endif
