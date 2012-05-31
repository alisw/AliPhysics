*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:50:55  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:11  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/07/11 17:14:49  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
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
