*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION MBYTET (MZ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    MBYTET          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

#if defined(CERNLIB_BUGLRSHFT)
      MBYTET = AND (IZW,
     +               NOT(ishft (LSHIFT(NOT(MZ),32-NZB), -33+IZP+NZB)))
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MBYTET = AND (IZW,
     +               NOT(lrshft (LSHIFT(NOT(MZ),32-NZB), 33-IZP-NZB)))
#endif
      RETURN
      END
