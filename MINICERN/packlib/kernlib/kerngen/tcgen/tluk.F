*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLUK (A,IASEP,NR,SIG,BETA)
C
C CERN PROGLIB# E230    TLUK            .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLUK (A,IASEP,NR,SIG,BETA)
C.
C.       COMPUTE TRANSFORMATION QUANTITIES.
C.       TLUK HAS BEEN MODIFIED FOR KINEMATICS.
C.
C.-------------------------------------------------------------------
C
      COMMON /SLATE/ DUM(37),I,JA,LL
      DIMENSION A(*)
C--
C
C--      COMPUTE MODULUS OF A GIVEN ROW IN A MATRIX AND FIND LAST
C--      NON-ZERO ELEMENT IN THAT ROW.
C
      SIG= 0.
      JA = 1
      LL = 0
C
      DO           10        I=1,NR
      IF     (A(JA).EQ.0.)             GO TO     10
      LL = I
      SIG= SIG + A(JA)* A(JA)
   10 JA = JA + IASEP
C
C--      FOR A ZERO ROW RETURN.
C
      NR = LL
      IF     (NR.EQ.0)                 RETURN
C--
C--      OTHERWISE TAKE THE MODULUS WITH SIGN OF FIRST ELEMENT OF THAT
C--      ROW. REDEFINE THAT FIRST ARGUMENT AND VALUE OF BETA.
C
      SIG  = SIGN (SQRT (SIG),A(1))
      BETA = A(1) + SIG
      A(1) = BETA
      BETA = 1. / (SIG * BETA)
      RETURN
      END
