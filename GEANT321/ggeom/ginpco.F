*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINPCO(X,P,IYES)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO PERFORM GINME FUNCTION FOR A POLYCONE.           *
C.    *        X  -  A VECTOR DIMENSION 3 GIVING THE POINT.            *
C.    *        P  -  THE VECTOR OF PARAMETERS, (1) = PHIMIN,           *
C.    *                   (2) = DPHI, (3) = NZ, (4) = Z(1),            *
C.    *                   (5) = RMIN(1), (6) = RMAX(1), (7) = Z(2),    *
C.    *                   (8) = RMIN(2), (9) = RMAX(2), ETC.           *
C.    *     IYES  -  RETURNED AS 1 IF POINT IS INSIDE, OTHERWISE 0.    *
C.    *                                                                *
C.    *    ==>Called by : GINME                                        *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcpoly.inc"
      DIMENSION P(9),X(3)
C.
C.                ---------------------------------
C.
      IYES=0
C
C                 FIRST CHECK TOTAL Z RANGE.
C
      NZ=P(3)
      IF(X(3).LT.P(4).OR.X(3).GT.P(NZ*3+1)) GO TO 999
C
C                 NOW CHECK R(SQUARED). FIRST COMPUTE IT.
C
      R2=X(1)*X(1)+X(2)*X(2)
C
C                 NOW FIND WHICH Z SECTION THE POINT IS IN.
C
      IZL=1
      IZH=NZ
      IZT=(NZ+1)/2
      IPL=IZL*3+1
      IPH=IZH*3+1
      IPT=IZT*3+1
C
   10 CONTINUE
      IF(IZH-IZL.LE.1) GO TO 40
C
      IF(X(3).GT.P(IPT)) GO TO 20
C
      IZH=IZT
      IPH=IPT
C
      GO TO 30
C
   20 CONTINUE
C
      IZL=IZT
      IPL=IPT
C
   30 CONTINUE
C
      IZT=(IZL+IZH)/2
      IPT=IZT*3+1
C
      GO TO 10
C
   40 CONTINUE
C
C          THE POINT IS IN THE SECTION BOUNDED BY THE IZL TH AND THE
C          IZH TH Z VALUES.
C
C          NOW COMPUTE RMIN AND RMAX AND TEST THE VALUE OF R SQUARED.
C
      IZSEC=IZH
      IF (P(IPL+2).EQ.0..AND.P(IPH+2).EQ.0.) GO TO 999
C
      IF ( (P(IPH).EQ.P(IPL)) .AND. (X(3).EQ.P(IPL)) ) THEN
         RMN = MAX( P(IPL+1), P(IPH+1) )
         RMX = MIN( P(IPL+2), P(IPH+2) )
      ELSE
         DZ=P(IPH)-P(IPL)
         DZ1=X(3)-P(IPL)
         RMN=(P(IPL+1)*(DZ-DZ1)+P(IPH+1)*DZ1)/DZ
         RMX=(P(IPL+2)*(DZ-DZ1)+P(IPH+2)*DZ1)/DZ
      END IF
C
      IF(R2.LT.RMN*RMN.OR.R2.GT.RMX*RMX) GO TO 999
C
C                 NOW CHECK PHI.
C
      PHI=P(1)+P(2)*0.5
      IF(X(2).NE.0.0.OR.X(1).NE.0) PHI=ATAN2(X(2),X(1))*RADDEG
      IF(PHI.LT.P(1)) PHI=PHI+360.0
      IF(PHI.LT.P(1).OR.PHI-P(1).GT.P(2)) GO TO 999
C
C           POINT IS IN POLYCONE.
C
      IYES=1
C
  999 CONTINUE
      END
