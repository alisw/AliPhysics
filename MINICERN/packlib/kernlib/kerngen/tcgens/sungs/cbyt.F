*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      SUBROUTINE CBYT (MZ,IZQ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    CBYT            .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
C
#if defined(CERNLIB_BUGLRSHFT)
      IZW = OR ( AND (IZW,
     +                    NOT(LSHIFT (ishft(NOT(0),NZB-32), IZP-1)))
     +,         ishft (LSHIFT(ishft(MZ,1-IZQ),32-NZB),-33+IZP+NZB) )
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      IZW = OR ( AND (IZW,
     +                    NOT(LSHIFT (lrshft(NOT(0),32-NZB), IZP-1)))
     +,         lrshft (LSHIFT(lrshft(MZ,IZQ-1),32-NZB), 33-IZP-NZB) )
#endif
      RETURN
      END
