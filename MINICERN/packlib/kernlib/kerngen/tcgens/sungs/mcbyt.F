*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION MCBYT (MZ,IZQ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    MCBYT           .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

#if defined(CERNLIB_BUGLRSHFT)
      MASK  = LSHIFT (ishft(-1,NZB-32), IZP-1)
      MCBYT = OR ( AND (IZW, NOT(MASK))
     +,            AND (MASK, LSHIFT(ishft(MZ,1-IZQ), IZP-1)))
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MASK  = LSHIFT (lrshft(-1,32-NZB), IZP-1)
      MCBYT = OR ( AND (IZW, NOT(MASK))
     +,            AND (MASK, LSHIFT(lrshft(MZ,IZQ-1), IZP-1)))
#endif
      RETURN
      END
