*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:35:40  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:30  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*
* q_sbit.inc
*
#if defined(CERNLIB_QMSUN)||defined(CERNLIB_QMLNX)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHIFT(1,IZP-1)) )
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHIFT(1,IZP-1)  )
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHIFT(1,IZP-1)) )
     +                             , LSHIFT(AND(MZ,1),IZP-1) )
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C)
      MSBIT0 (IZW,IZP)   = AND (IZW, NOT(LSHIFT(1,IZP-1)) )
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHIFT(1,IZP-1)  )
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHIFT(1,IZP-1)) )
     +                              , LSHIFT(AND(MZ,1),IZP-1))
#elif defined(CERNLIB_QMOS9)||defined(CERNLIB_QMHPX)||(defined(CERNLIB_QMDOS)&&(!defined(CERNLIB_QF_F2C))) || defined(CERNLIB_QFMSOFT)
      MSBIT0 (IZW,IZP)   = IAND (IZW, NOT(ISHFT(1,IZP-1)) )
      MSBIT1 (IZW,IZP)   = IOR  (IZW,     ISHFT(1,IZP-1)  )
      MSBIT (MZ,IZW,IZP) = IOR (IAND (IZW, NOT(ISHFT(1,IZP-1)) )
     +                              , ISHFT(IAND(MZ,1),IZP-1))
#endif
