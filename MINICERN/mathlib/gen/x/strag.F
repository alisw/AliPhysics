*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION STRAG(PINV,EM)
      DIMENSION ENT(15),TAB(15)
      DATA  ENT/.0,.001,.00189,.00355,.00675,.0128,.0282,.0460,.0875,.16
     +   5,.312,.593,1.012,5.,100000./
      DATA  TAB/.1,.01573,.01445,.0134,.0124,.0117,.0112,.0107,.0102,.00
     +   97,.0092,.0087,.0082,0.,0./
      BI     =SQRT(1.+1./(PINV*EM)**2)-1.
      DO 1 I=2,15
      IF(ENT(I).GT.BI)    GO TO 2
    1 CONTINUE
    2 A=(BI-ENT(I-1))*(TAB(I)-TAB(I-1))/(ENT(I)-ENT(I-1))+TAB(I-1)
      STRAG=SQRT(.938213/EM)*A
      RETURN
      END
