*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SPLAN0(B0,X,X1,X2,ND)
#include "gen/imp64.inc"
      DIMENSION B0(1)

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   TEST-ROUTINE  SPLAN0
*   FOR COMPUTING  1-DIMENSIONAL  B-SPLINES OF DEGREE  K = 0
*   BY MEANS OF ANALYTIC FORMULAE
*
************************************************************************

      B0(1)=0
      IF(ND .GT. 0) RETURN
      IF(ND .EQ. 0) THEN
       IF(X1 .LE. X .AND. X .LE. X2) B0(1)=1
       RETURN
      ENDIF
      IF(X1 .LT. X) B0(1)=MIN(1D0,X-X1)
      RETURN
      END



