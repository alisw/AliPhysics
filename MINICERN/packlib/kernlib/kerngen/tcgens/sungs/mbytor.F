*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION MBYTOR (MZ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    MBYTOR          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

#if defined(CERNLIB_BUGLRSHFT)
      MBYTOR = OR (IZW, ishft (LSHIFT(MZ,32-NZB), -33+IZP+NZB) )
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MBYTOR = OR (IZW, lrshft (LSHIFT(MZ,32-NZB), 33-IZP-NZB) )
#endif
      RETURN
      END
