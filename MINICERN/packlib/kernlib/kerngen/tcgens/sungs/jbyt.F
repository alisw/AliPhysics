*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION JBYT (IZW,IZP,NZB)
C
C CERN PROGLIB# M421    JBYT            .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
#if defined(CERNLIB_BUGLRSHFT)
      JBYT = ishft (LSHIFT(IZW,33-IZP-NZB),-32+NZB)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      JBYT = lrshft (LSHIFT(IZW,33-IZP-NZB), 32-NZB)
#endif
      RETURN
      END
