*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION JRSBYT (IT,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    JRSBYT          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

#if defined(CERNLIB_BUGLRSHFT)
      MASK   = LSHIFT (ishft(-1,-32+NZB), IZP-1)
      JRSBYT = ishft (AND(MASK,IZW), 1-IZP)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MASK   = LSHIFT (lrshft(-1,32-NZB), IZP-1)
      JRSBYT = lrshft (AND(MASK,IZW), IZP-1)
#endif
      IZW    = OR (AND(NOT(MASK),IZW), AND(MASK,LSHIFT(IT,IZP-1)))
      RETURN
      END
