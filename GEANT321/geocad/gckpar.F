*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GCKPAR (ISH,NPAR,PAR)
*
#include "geant321/gconsp.inc"
*
      DIMENSION PAR(50)
*
      IF (ISH.EQ.4) THEN
* ***    trap
         NPAR = 11
         PH = 0.
         IF (PAR(2).NE.0.) PH = ATAN2(PAR(3),PAR(2))*RADDEG
         TT = SQRT(PAR(2)**2+PAR(3)**2)
         PAR(2) = ATAN(TT)*RADDEG
         IF (PH.LT.0.0) PH = PH + 360.0
         PAR(3) = PH
         PAR(7) = ATAN(PAR(7))*RADDEG
         IF (PAR(7).GT.90.0) PAR(7) = PAR(7)-180.0
         PAR(11)= ATAN(PAR(11))*RADDEG
         IF (PAR(11).GT.90.0) PAR(11) = PAR(11)-180.0
*
      ELSE IF (ISH.EQ.10) THEN
* ***    para
         PH = 0.
         IF (PAR(5).NE.0.) PH = ATAN2(PAR(6),PAR(5))*RADDEG
         TT = SQRT(PAR(5)**2+PAR(6)**2)
         PAR(4) = ATAN(PAR(4))*RADDEG
         IF (PAR(4).GT.90.0) PAR(4) = PAR(4)-180.0
         PAR(5) = ATAN(TT)*RADDEG
         IF (PH.LT.0.0) PH = PH + 360.0
         PAR(6) = PH
*
      ELSE IF (ISH.EQ.28) THEN
* ***    gtra
         NPAR = 12
      ENDIF
*
      END
