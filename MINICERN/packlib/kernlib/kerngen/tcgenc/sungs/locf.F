*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:38  mclareni
* Kernlib
*
*
      FUNCTION LOCF (AM)
C
C CERN PROGLIB# N100    LOCF            .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
      DIMENSION AM(9)

      LLAD = LOC(AM(1))
#if defined(CERNLIB_BUGLRSHFT)
      LOCF = ishft (LLAD, -2)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      LOCF = lrshft (LLAD, 2)
#endif
      RETURN
      END
