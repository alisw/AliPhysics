*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:45  mclareni
* Mathlib gen
*
*
*FCA :          Fri Mar 26 17:27:50 CET 1999 by  Federico Carminati
*               removed VZERO
#include "gen/pilot.h"
      SUBROUTINE U100FC (I,K)
      DIMENSION K(40), L(40)
      DO KZERO=1,40
         K(KZERO)=0
      ENDDO
*      CALL VZERO (K, 40)
      IF(I)  4,3,2
    4 K(1)=-77
      RETURN
    2 DO 5  JK=1,I
      CALL U100PR(JK,L,M)
      DO 5 JM=1,40
    5 K(JM)=K(JM)+L(JM)
    3 RETURN
      END
