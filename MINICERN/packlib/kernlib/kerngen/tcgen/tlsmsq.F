*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLSMSQ (B,L,M,F)
C
C CERN PROGLIB# E230    TLSMSQ          .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLSMSQ (B,L,M,F)
C.
C.       CALCULATES THE SUM OF SQUARES OF RESIDUALS.
C.
C.-------------------------------------------------------------------
C
      COMMON /SLATE/ DUM(38),I,JB
      DIMENSION      B(*)
C
      F  = 0.
      JB = 1
      DO           10        I=1,M
      F  = F + B(JB)*B(JB)
   10 JB = JB + L
C
      RETURN
      END
