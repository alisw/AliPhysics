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
* Revision 1.1.1.1  1996/02/15 17:49:19  mclareni
* Kernlib
*
*
*
* q_cbyt.inc
*
#if defined(CERNLIB_QMLNX)
      MCBYT (MZ,IZQ,IZW,IZP,NZB) = OR (
     +       AND (IZW, NOT(LSHIFT (ISHFTR(NOT(0),32-NZB),IZP-1)))
     +      ,ISHFTR (LSHIFT(LSHIFT(MZ,IZQ-1),32-NZB), 33-IZP-NZB))
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C)
C      MCBYT (MZ,IZQ,IZW,IZP,NZB) = OR (
C     +          AND (IZW, NOT(LSHIFT (ishftr(NOT(0),32-NZB), IZP-1)))
C     +        , ishftr (LSHIFT(LSHIFT(MZ,-IZQ+1),32-NZB), 33-IZP-NZB) )
#elif defined(CERNLIB_QMOS9)||defined(CERNLIB_QMHPX)||(defined(CERNLIB_QMDOS)&&(!defined(CERNLIB_QF_F2C)))
      MCBYT (MZ,IZQ,IZW,IZP,NZB) = IOR (
     +          IAND (IZW, NOT(ISHFT (ISHFT(NOT(0),-32+NZB), IZP-1)))
     +        , ISHFT (ISHFT(ISHFT(MZ,-IZQ+1),32-NZB), -33+IZP+NZB) )
#endif
