*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:35:38  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:28  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:19  mclareni
* Kernlib
*
*
*    IAND / IOR / IEOR / NOT  standard?
*
* q_andor.inc
*
#if defined(CERNLIB_QMSUN)||defined(CERNLIB_QMLNX)
      IAND (IZV,IZW) = AND (IZV, IZW)
      IOR  (IZV,IZW) =  OR (IZV, IZW)
      IEOR (IZV,IZW) = XOR (IZV, IZW)
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF2C)
      IAND(IZV,IZW)  = AND(IZV,IZW)
       IOR(IZV,IZW)  =  OR(IZV,IZW)
      IEOR(IZV,IZW)  = XOR(IZV,IZW)
#else
#endif
