*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      SUBROUTINE SBYT (MZ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    SBYT            .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

#if defined(CERNLIB_BUGLRSHFT)
      MASK = LSHIFT (ishft(-1,NZB-32), IZP-1)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MASK = LSHIFT (lrshft(-1,32-NZB), IZP-1)
#endif
      IZW  = OR (AND(NOT(MASK),IZW), AND(MASK,LSHIFT(MZ,IZP-1)))
      RETURN
      END
