*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*
* q_isabit.inc
*
#if defined(CERNLIB_QMLNX)||(defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C))
      BTEST(IZW,IZP) = AND (IZW,LSHIFT(1,IZP)) .NE. 0
      IBSET(IZW,IZP) =  OR (IZW,LSHIFT(1,IZP))
      IBCLR(IZW,IZP) = AND (IZW,NOT(LSHIFT(1,IZP)))
      IBITS (IZW,IZP,NZB) = ISHFTR (LSHIFT(IZW,32-IZP-NZB), 32-NZB)
#endif
