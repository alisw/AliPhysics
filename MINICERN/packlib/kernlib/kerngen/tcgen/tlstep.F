*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLSTEP (A,B,IASEP,IBSEP,NR,NC,BETA)
C
C CERN PROGLIB# E230    TLSTEP          .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLSTEP (A,B,IASEP,IBSEP,NR,NC,BETA)
C.
C.       PERFORMS AN ORTHOGONAL TRANSFORMATION.
C.
C.-------------------------------------------------------------------
C
      COMMON /SLATE/ DUM(34),H,I,IB,J,JA,JB
      DIMENSION      A(*), B(*)
C
      IB = 0
      DO           30        J=1,NC
      IB = IB + 1
      H  = 0.
      JA = 1
      JB = IB
C
      DO           10        I=1,NR
      H  = H + A(JA) * B(JB)
      JA = JA +IASEP
   10 JB = JB + IBSEP
C
      H  = H * BETA
      JA = 1
      JB = IB
      DO           20        I=1,NR
      B(JB) = B(JB) - A(JA) * H
      JA = JA +IASEP
   20 JB = JB + IBSEP
   30 CONTINUE
C
      RETURN
      END
