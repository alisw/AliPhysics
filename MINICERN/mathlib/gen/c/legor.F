*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LEGOR
#include "legbl.inc"
     1,RGAM
      A=0.5*(1.0+VV)
      B=0.0
      C=0.5
      GR=RGAM(A,B,C)
      C=U*PI*VV*0.5
      A=EXP(C)
      B=EXP(-C)
      IF(NFRIG)9,10,11
    9 QQ=PISR*0.5*GR*U*A
      GOTO 12
   10 QQ=U*PISR/4.0*GR*(A-B)
      GOTO 12
   11 QQ=-PISR/2.0*GR*B*U
   12 PP=GR*(A+B)/(2.0*PISR)
      RETURN
      END
