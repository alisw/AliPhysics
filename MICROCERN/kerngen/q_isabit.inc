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
