*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:21  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ARITHM (Y)
C-----ARITHM SETS Y AND CALCULATES DY/DA
#include "d510si.inc"
#include "d510ui.inc"
C-----10.*MAXIMUM RELATIVE PRECISION ON CDC 6000
      DATA RP/1.E-14/
      Y=FUNCT(X)
      DO 6 I=1,NA
      DF(I)=0.
      IF (PL0(I)) 6,6,1
 1    AI=A(I)
      HI=0.01*PL0(I)
      PI=RP*ABS(AI)
      IF (HI.LE.PI) HI=PI
      A(I)=AI+HI
      IF (A(I)-AMX(I)) 5,5,2
 2    A(I)=AI-HI
      HI=-HI
      IF (A(I)-AMN(I)) 3,5,5
 3    A(I)=AMX(I)
      HI=AMX(I)-AI
      IF (AMN(I)-AI+HI) 4,5,5
 4    A(I)=AMN(I)
      HI=AMN(I)-AI
 5    DF(I)=(FUNCT(X)-Y)/HI
      A(I)=AI
 6    CONTINUE
      RETURN
      END
