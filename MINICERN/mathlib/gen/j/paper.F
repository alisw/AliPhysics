*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE PAPER(IMAGEW,IMAGEH)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PURPOSE                                                            C
C     TO DEFINE HOW LARGE PICTURES SHALL BE PRODUCED BY SUBSEQUENT     C
C     CALLS TO THE SUBROUTINE CONT.                                    C
C   USAGE                                                              C
C     CALL PAPER(IMAGEW,IMAGEH)                                        C
C   PARAMETERS                                                         C
C     IMAGEW - THE IMAGE WIDTH IN NUMBER OF LETTERS ON A LINE          C
C     IMAGEH - THE IMAGE HEIGT IN NUMBER OF LINES                      C
C   REMARKS                                                            C
C     IMAGEW AND IMAGEH SHOULD NOT EXEED CERTAIN MACHINE-DEPENDANT     C
C     LIMITS,  USUALLY 120 AND 62, RESPECTIVELY.
C
      COMMON /J509C1/XBL,YBL,XUR,YUR,NX(6),NY(6),LETT1,LETT2,LINE1,LINE2
C
      CALL J509BD
      LETT2 = IMAGEW - 1
      LINE2 = IMAGEH - 3
      IF (LETT2 .GT. 119)  LETT2 = 119
      RETURN
      END
