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
      SUBROUTINE GINME(X,SH,P,IYES)
C.
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     SUBROUTINE TO COMPUTE WHETHER THE POINT X IS IN THE SHAPE  *
C.    *     NUMBER SH WITH PARAMETERS P.                               *
C.    *     IYES IS SET TO 1 IF SO; TO ZERO IF NOT.                    *
C.    *                                                                *
C.    *    ==>Called by : GFTRAC, GINVOL, GMEDIA, GMEPOS, GSTRAC,      *
C.    *                   GTMEDI                                       *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      DIMENSION X(3),P(*)
C.
C.    ------------------------------------------------------------------
C.
      IYES=0
      ISH=SH
*
      IF(ISH.EQ.5) THEN
*
*          Tube axis = Z axis
         IF (ABS(X(3)).GT.P(3)) GO TO 999
         R2=X(1)**2+X(2)**2
         IF(R2.LT.P(1)**2) GO TO 999
         IF(R2.GT.P(2)**2) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.6) THEN
*
*          Tube segment
         IF (ABS(X(3)).GT.P(3)) GO TO 999
         R2=X(1)**2+X(2)**2
         IF(R2.LT.P(1)**2) GO TO 999
         IF(R2.GT.P(2)**2) GO TO 999
         IYES=1
         IF (R2.LT.1.0E - 6) GO TO 999
         IYES=0
         PHI = ATAN2(X(2),X(1)) * RADDEG
         DP=P(5)-P(4)
         IF(DP.LE.0.0) DP=DP+360.0
         DDP=PHI-P(4)
         SG = SIGN(1.0,DDP)
         DDP = MOD( ABS(DDP), 360. )
         IF(SG.LE.0.) DDP = 360.-DDP
         IF(DDP.GT.DP) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.1) THEN
*
*            Rectilinear box
         IF (ABS(X(1)).GT.P(1)) GO TO 999
         IF (ABS(X(2)).GT.P(2)) GO TO 999
         IF (ABS(X(3)).GT.P(3)) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.2) THEN
*
*          Rectilinear trapezoidal section with only the X
*          thickness changing with z
         DZ =P(4)
         IF (ABS(X(3)).GT.DZ) GO TO 999
         DY = P(3)
         IF (ABS(X(2)).GT.DY) GO TO 999
         DX =0.5*(P(2)*(X(3)+DZ)+P(1)*(DZ-X(3)))/DZ
         IF (ABS(X(1)).GT.DX) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.3) THEN
*
*          Rectilinear trapezoidal section with both X and Y
*          thicknesses varying with Z
         DZ =P(5)
         IF (ABS(X(3)).GT.DZ) GO TO 999
         DX =0.5*(P(2)*(X(3)+DZ)+P(1)*(DZ-X(3)))/DZ
         IF (ABS(X(1)).GT.DX) GO TO 999
         DY = 0.5*(P(4)*(X(3)+DZ)+P(3)*(DZ-X(3)))/DZ
         IF (ABS(X(2)).GT.DY) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.9) THEN
*
*          Spherical segment
         R2 = X(1)**2 + X(2)**2 + X(3)**2
         IF (R2.LT.P(1)**2.OR.R2.GT.P(2)**2) GO TO 999
         IYES=1
         IF (R2.LT.1.0E -12) GO TO 999
         IYES=0
         IF (X(1).EQ.0..AND.X(2).EQ.0.) THEN
            PHI = 0.
         ELSE
            PHI = ATAN2(X(2),X(1))*RADDEG
         ENDIF
         DP=P(6)-P(5)
         IF(DP.LE.0.0) DP=DP+360.0
         DDP=PHI-P(5)
         SG = SIGN(1.0,DDP)
         DDP = MOD( ABS(DDP), 360. )
         IF(SG.LE.0.) DDP = 360.-DDP
         IF(DDP.GT.DP) GO TO 999
         R2 = X(1)**2 + X(2)**2
         IF(R2.GT.0.0) R2=SQRT(R2)
         THETA = ATAN2(R2,X(3)) * RADDEG
         IF (THETA.LT.P(3).OR.THETA.GT.P(4)) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.4) THEN
*
*          General trapezoidal section
         CALL GINTRP(X,P,IYES)
*
      ELSEIF(ISH.EQ.7) THEN
*
*          Conical tube
         DZ =P(1)
         IF (ABS(X(3)).GT.DZ) GO TO 999
         R2 =X(1)**2 + X(2)**2
         RL = 0.5*(P(4)*(X(3)+DZ)+P(2)*(DZ-X(3)))/DZ
         RH = 0.5*(P(5)*(X(3)+DZ)+P(3)*(DZ-X(3)))/DZ
         IF (R2.LT.RL**2.OR.R2.GT.RH**2) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.8) THEN
*
*          Conical tube segment
         DZ =P(1)
         IF (ABS(X(3)).GT.DZ) GO TO 999
         R2 =X(1)**2 + X(2)**2
         RL = 0.5*(P(4)*(X(3)+DZ)+P(2)*(DZ-X(3)))/DZ
         RH = 0.5*(P(5)*(X(3)+DZ)+P(3)*(DZ-X(3)))/DZ
         IF (R2.LT.RL**2.OR.R2.GT.RH**2) GO TO 999
         IYES=1
         IF (R2.LT.1.0E - 6) GO TO 999
         IYES=0
         PHI = ATAN2(X(2),X(1)) * RADDEG
         DP=P(7)-P(6)
         IF(DP.LE.0.0) DP=DP+360.0
         DDP=PHI-P(6)
         SG = SIGN(1.0,DDP)
         DDP = MOD( ABS(DDP), 360. )
         IF(SG.LE.0.) DDP = 360.-DDP
         IF(DDP.GT.DP) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.10) THEN
*
*             Parallelepiped
         CALL GINPAR(X,P,IYES)
*
      ELSEIF(ISH.EQ.11) THEN
*
*            Polygon
         CALL GINPGO(X,P,IYES)
*
      ELSEIF(ISH.EQ.12) THEN
*
*            Polycone
         CALL GINPCO(X,P,IYES)
*
      ELSEIF (ISH.EQ.13) THEN
*
*           Elliptical tube
         IF (ABS(X(3)).GT.P(3)) GOTO 999
         A2=P(1)**2
         B2=P(2)**2
         R2=X(1)**2/A2+X(2)**2/B2
         IF (R2.GT.1.0) GOTO 999
         IYES=1
*
      ELSEIF(ISH.EQ.14) THEN
*
*       HYPErboloid axis = Z axis.
         IF (ABS(X(3)) .GT. P(3)) GO TO 999
         R2 = X(1)**2 + X(2)**2
         DR2 = (X(3) * TAN(P(4)*DEGRAD))**2
         IF (R2 .LT. P(1)**2 + DR2) GO TO 999
         IF (R2 .GT. P(2)**2 + DR2) GO TO 999
         IYES=1
*
      ELSEIF(ISH.EQ.28) THEN
*
*            General twisted trapezoid.
         CALL GINGTR(X,P,IYES)
*
      ELSEIF( ISH.EQ.NSCTUB) THEN
*
*            Cut tube.
         CALL GINCTU(X,P,IYES)
*
      ELSE
*
*            User shape ?
#if !defined(CERNLIB_USRJMP)
         CALL GUINME(X,SH,P,IYES)
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT4(JUINME,X,SH,P,IYES)
#endif
*
         IF(IYES.LT.0) THEN
*
*         Shape number not supported by system nor by user
            WRITE(CHMAIL,10000) ISH
            CALL GMAIL(0,0)
         ENDIF
      ENDIF
*
10000 FORMAT(' *** GINME *** Shape number ',I5,' not yet implemented')
  999 END
