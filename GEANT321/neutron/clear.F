*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE CLEAR(L,L1,L2)
C       THIS ROUTINE ZEROS ARRAY L FROM
C       STARTING POINT L1 TO ENDING POINT L2
      DIMENSION L(*)
      IF(L2-L1.LT.0)GO TO 20
      DO 10 I=L1,L2
   10 L(I)=0
   20 RETURN
      END
