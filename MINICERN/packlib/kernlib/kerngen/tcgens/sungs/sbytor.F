*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      SUBROUTINE SBYTOR (MZ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    SBYTOR          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

#if defined(CERNLIB_BUGLRSHFT)
      IZW = OR (IZW, LSHIFT (AND(ishft(-1,NZB-32),MZ), IZP-1))
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      IZW = OR (IZW, LSHIFT (AND(lrshft(-1,32-NZB),MZ), IZP-1))
#endif
      RETURN
      END
