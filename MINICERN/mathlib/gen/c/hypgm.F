*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE HYPGM(A,B,C,Z,H,ACC,NCVG)
      COMPLEX A,B,C,Z,H,AA,BB,CC,ZZ,ADD,DD,HH
      NCVG=0
      ZZ=Z
      HH=(1.0,0.0)
      AA=A
      BB=B
      CC=C
      ADD=(1.0,0.0)
      DD=(1.0,0.0)
      DO 21 I=1,50
      ADD=ADD*ZZ*AA/DD*BB/CC
      HH=HH+ADD
      IF(CRIT(HH,ADD,ACC))3,4,4
    4 AA=AA+1.0
      BB=BB+1.0
      CC=CC+1.0
      DD=DD+1.0
   21 CONTINUE
      NCVG=1
    3 H=HH
      RETURN
      END
