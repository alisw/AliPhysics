*
* $Id$
*
* $Log$
* Revision 1.3  1997/09/02 14:26:58  mclareni
* WINNT correction
*
* Revision 1.2  1997/02/04 17:35:40  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:31  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*
* q_sbyt.inc
*
#if defined(CERNLIB_QMLNX)
      MSBYT (MZ,IZW,IZP,NZB) = OR (
     +       AND (IZW, NOT(LSHIFT (ISHFTR(NOT(0),32-NZB),IZP-1)))
     +      ,ISHFTR (LSHIFT(MZ,32-NZB), 33-IZP-NZB) )
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C)
C      MSBYT (MZ,IZW,IZP,NZB) = OR (
C     +          AND (IZW, NOT(LSHIFT (ISHFTR(NOT(0),32-NZB), IZP-1)))
C     +        , ISHFTR (LSHIFT(MZ,32-NZB), 33-IZP-NZB) )
#elif defined(CERNLIB_QMOS9)||defined(CERNLIB_QMHPX)||(defined(CERNLIB_QMDOS)&&(!defined(CERNLIB_QF_F2C))) || defined(CERNLIB_QFMSOFT) || defined(CERNLIB_QFDEC)
      MSBYT (MZ,IZW,IZP,NZB) = IOR (
     +          IAND (IZW, NOT(ISHFT (ISHFT(NOT(0),-32+NZB), IZP-1)))
     +        , ISHFT (ISHFT(MZ,32-NZB), -33+IZP+NZB) )
#endif
