*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNVAX  2.31  890323
C ORIG. 16/03/89  JZ
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF
C
      COMMON /SLATE/ MCH(40)

      CALL UBLOW (
     + 40HKERNVAX  2.47  1995/07/18 13.58     CERN
     +,     MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
