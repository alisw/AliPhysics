*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFLPHI(ISH,IROT,DX,PARS,CL,CH,IERR)
C.
C.    **********************************************************
C.    *                                                        *
C.    *    ROUTINE TO FIND THE PHI LIMITS OF THE OBJECT SHAPE  *
C.    *    ISH ROTATED BY MATRIX IROT AND DISPLACED BY VECTOR  *
C.    *    DX. IT HAS NPAR PARAMTERS IN THE ARRAY PARS. THE    *
C.    *    LOWER LIMIT IS RETURNED IN CL AND THE HIGHER IN CH. *
C.    *    NOTE THE OBJECT IS CONTAINED IN THE RANGE OF        *
C.    *    INCREASING PHI FROM CL TO CH THOUGH CL AND CH ARE   *
C.    *    FORCED TO LIE IN THE RANGE 0.0 TO 360.0 SO THAT THE *
C.    *    VALUE OF CL CAN BE HIGHER THAN THAT OF CH.          *
C.    *                                                        *
C.    *    ==>Called by : GFCLIM                               *
C.    *         Author  A.McPherson  *********                 *
C.    *                                                        *
C.    **********************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
      DIMENSION DX(3),PARS(11),X(3),XT(3)
C.
C.          -------------------------------------------
C.
      IERR=1
C
      DXS=DX(1)*DX(1)+DX(2)*DX(2)
      IF(DXS.GT.0.0) DXS=SQRT(DXS)
C
      IF(ISH.GT.4.AND.ISH.NE.10.AND.ISH.NE.28) GO TO 40
C
C
C           CUBOIDS, TRAPEZOIDS, PARALLELEPIPEDS.
C
      IERR=0
      CL=0.0
      CH=360.0
C
C           IF IN DOUBT SET IT TO FULL RANGE.
C
      IF(DXS.LE.0.0) GO TO 999
C
      PHC=90.
      IF(DX(1).NE.0.)PHC=ATAN2(DX(2),DX(1))*RADDEG
      IF(PHC.LT.0.0) PHC=PHC+360
      PL=0.0
      PH=0.0
C
      DO 30 IP=1,8
C
C           THIS IS A LOOP OVER THE 8 CORNERS.
C           FIRST FIND THE LOCAL COORDINATES.
C
      IF(ISH.EQ.28) THEN
C
C            General twisted trapezoid.
C
         IL=(IP+1)/2
         I0=IL*4+11
         IS=(IP-IL*2)*2+1
         X(3)=PARS(1)*IS
         X(1)=PARS(I0)+PARS(I0+2)*X(3)
         X(2)=PARS(I0+1)+PARS(I0+3)*X(3)
         GO TO 20
C
      ENDIF
C
      IP3=ISH+2
      IF(ISH.EQ.10) IP3=3
      IF(ISH.EQ.4) IP3=1
      X(3)=PARS(IP3)
      IF(IP.LE.4) X(3)=-X(3)
      IP2=3
      IF(ISH.GT.2.AND.X(3).GT.0.0) IP2=4
      IF(ISH.EQ.1.OR.ISH.EQ.10) IP2=2
      IF(ISH.EQ.4) IP2=4
      IF(ISH.EQ.4.AND.X(3).GT.0.0) IP2=8
      X(2)=PARS(IP2)
      IF(MOD(IP+3,4).LT.2) X(2)=-X(2)
      IP1=1
      IF(ISH.NE.1.AND.ISH.NE.10.AND.X(3).GT.0.0) IP1=2
      IF(ISH.EQ.4) IP1=5
      IF(ISH.EQ.4.AND.X(3).GT.0.0) IP1=IP1+4
      IF(ISH.EQ.4.AND.X(2).GT.0.0) IP1=IP1+1
      X(1)=PARS(IP1)
      IF(MOD(IP,2).EQ.1) X(1)=-X(1)
C
      IF(ISH.NE.10) GO TO 10
      X(1)=X(1)+X(2)*PARS(4)+X(3)*PARS(5)
      X(2)=X(2)+X(3)*PARS(6)
   10 CONTINUE
C
      IF(ISH.NE.4) GO TO 20
      IP4=7
      IF(X(3).GT.0.0) IP4=11
      X(1)=X(1)+X(2)*PARS(IP4)+X(3)*PARS(2)
      X(2)=X(2)+X(3)*PARS(3)
   20 CONTINUE
C
C          ROTATE.
C
      JROT=LQ(JROTM-IROT)
      XT(1)=X(1)
      XT(2)=X(2)
      XT(3)=X(3)
      IF(IROT.NE.0) CALL GINROT(X,Q(JROT+1),XT)
C
      XPT=DXS+(DX(1)*XT(1)+DX(2)*XT(2))/DXS
      YPT=(DX(1)*XT(2)-DX(2)*XT(1))/DXS
C
      IF(YPT.EQ.0.0.AND.XPT.EQ.0.0) GO TO 999
      P=ATAN2(YPT,XPT)
      IF(P.GT.PI) P=P-PI*2.0
      IF(P.LT.PL) PL=P
      IF(P.GT.PH) PH=P
C
C
   30 CONTINUE
C
C
      IF(PH-PL.GT.PI) GO TO 999
      CL=PHC+PL*RADDEG
      CH=PHC+PH*RADDEG
C
      SG = SIGN(1.0,CL)
      CL = MOD( ABS(CL),360.0 )
      IF(SG.LE.0.0) CL=360.-CL
      SG=SIGN(1.0,CH)
      CH = MOD( ABS(CH),360.0 )
      IF(SG.LE.0.0) CH=360.-CH
C
      GO TO 999
C
   40 CONTINUE
      IF(ISH.GT.8.AND.ISH.NE.NSCTUB.AND.ISH.NE.13.AND.ISH.NE.14)GO TO 80
C
C             TUBES AND CONES.
C
      IERR=0
      CL=0.0
      CH=360.0
C
C             WHEN IN DOUBT SET TO FULL RANGE.
C
      RM=PARS(2)
      IF(ISH.LE.6.OR.ISH.EQ.NSCTUB) GO TO 50
**
      IF(ISH.EQ.13) THEN
**
**       approxime to a cylinder whit radius
**       equal to the ellipse major axis
**
         IF(PARS(1).GT.RM) RM=PARS(1)
         GOTO 50
      ENDIF
      IF(ISH.EQ.14) THEN
        RM = SQRT(PARS(2)**2+(PARS(3)*TAN(PARS(4)*DEGRAD))**2)
        GO TO 50
      ENDIF
      RM=PARS(3)
      IF(PARS(5).GT.RM) RM=PARS(5)
C
   50 CONTINUE
C
      IF(DXS.GT.RM) GO TO 70
      IF(ISH.EQ.5.OR.ISH.EQ.7.OR.ISH.EQ.14) GO TO 999
      IF(ISH.EQ.13) GOTO 999
*                 Here we treat the CONS
      IF(IROT.EQ.0) THEN
*                 This is the simple case, no rotation
*                 Compute the position of the limits on
*                 the X-Y plane.
      PHIMIN=PARS(4)*DEGRAD
      PHIMAX=PARS(5)*DEGRAD
      DDX1 = DX(1)+RM*COS(PHIMIN)
      DDY1 = DX(2)+RM*SIN(PHIMIN)
      DDX2 = DX(1)+RM*COS(PHIMAX)
      DDY2 = DX(2)+RM*SIN(PHIMAX)
      CL = ATAN2(DDY1,DDX1)*RADDEG
      CH = ATAN2(DDY2,DDX2)*RADDEG
      ELSE
*                 Rotated tubes might be more difficult
*                 Just leave it for later
      CONTINUE
      ENDIF
      IF(ISH.LE.6) GO TO 60
      CL=PARS(6)
      CH=PARS(7)
   60 CONTINUE
C
      IF(IROT.EQ.0) GO TO 65
      JROT=LQ(JROTM-IROT)
      IF(Q(JROT+15).NE.0.0.AND.Q(JROT+15).NE.180.0) GO TO 999
C
      PHX=Q(JROT+12)
      PHY=Q(JROT+14)
      IF(PHY.LT.PHX) PHY=PHY+360.0
      ISPH=1
      IF(PHY-PHX.GT.180.0) ISPH=-1
      CL=ISPH*CL+PHX
      CH=ISPH*CH+PHX
      IF(ISPH.EQ.1) GO TO 65
      CHT=CH
      CH=CL
      CL=CHT
C
   65 CONTINUE
C
      SG=SIGN(1.0,CL)
      CL = MOD( ABS(CL),360.0 )
      IF(SG.LE.0.0) CL=360.-CL
      SG=SIGN(1.0,CH)
      CH = MOD( ABS(CH),360.0 )
      IF(SG.LE.0.0) CH=360.-CH
C
      GO TO 999
C
   70 CONTINUE
C
C            DISPLACEMENT GREATER THAN MAXIMUM RADIUS SO
C            ASSUME COMPLETE TUBE TO GENERATE 'WORST CASE'.
C
      DZ=PARS(3)
      IF(ISH.EQ.NSCTUB) THEN
        S1 = (1.0-PARS(8))*(1.0+PARS(8))
        IF( S1 .GT. 0.0) S1 = SQRT(S1)
        S2 = (1.0-PARS(11))*(1.0+PARS(11))
        IF( S2 .GT. 0.0) S2 = SQRT(S2)
        IF( S2 .GT. S1 ) S1 = S2
        DZ = DZ+RM*S1
      ELSEIF(ISH.GT.6.AND.ISH.NE.13.AND.ISH.NE.14) THEN
        DZ=PARS(1)
      ENDIF
C
      X(1)=0.0
      X(2)=0.0
      X(3)=1.0
C
C                    LOCAL Z AXIS.
C
      JROT=LQ(JROTM-IROT)
      XT(1)=X(1)
      XT(2)=X(2)
      XT(3)=X(3)
      IF(IROT.NE.0) CALL GINROT(X,Q(JROT+1),XT)
C
      COST=ABS(DX(1)*XT(1)+DX(2)*XT(2))
      COS2=ABS(DX(1)*XT(2)-DX(2)*XT(1))
      SINT=(DXS+COST)*(DXS-COST)
      SIN2=(DXS+COS2)*(DXS-COS2)
      IF(SINT.GT.0.0) SINT=SQRT(SINT)
      IF(SIN2.GT.0.0) SIN2=SQRT(SIN2)
C
      XPT=DXS-(COST*DZ+SINT*RM)/DXS
C
      IF(XPT.LE.0.0) GO TO 999
      YPT=(SIN2*RM+COS2*DZ)/DXS
      DP=ATAN(YPT/XPT)
C
      P0=ATAN2(DX(2),DX(1))
      CL=(P0-DP)*RADDEG
      CH=(P0+DP)*RADDEG
C
      SG=SIGN(1.0,CL)
      CL = MOD( ABS(CL),360.0 )
      IF(SG.LE.0.0) CL=360.-CL
      SG=SIGN(1.0,CH)
      CH = MOD( ABS(CH),360.0 )
      IF(SG.LE.0.0) CH=360.-CH
C
      GO TO 999
C
   80 CONTINUE
      IF(ISH.GT.9) GO TO 999
C
C               SPHERE.
C
      IERR=0
      CL=0.0
      CH=360.0
C
      IF(IROT.NE.0.OR.DXS.GT.0.0) GO TO 90
C
C          UNROTATED AND CENTERED.
C
      CL=PARS(5)
      CH=PARS(6)
C
      SG=SIGN(1.0,CL)
      CL = MOD( ABS(CL),360.0 )
      IF(SG.LE.0.0) CL=360.-CL
      SG=SIGN(1.0,CH)
      CH = MOD( ABS(CH),360.0 )
      IF(SG.LE.0.0) CH=360.-CH
C
      GO TO 999
C
   90 CONTINUE
C
C            ROTATED OR NOT CENTERED.
C
      IF(DXS.LT.PARS(2)) GO TO 999
      P0=ATAN2(DX(2),DX(1))
      DP=ASIN(PARS(2)/DXS)
      CL=(P0-DP)*RADDEG
      CH=(P0+DP)*RADDEG
C
      SG=SIGN(1.0,CL)
      CL = MOD( ABS(CL),360.0 )
      IF(SG.LE.0.0) CL=360.-CL
      SG=SIGN(1.0,CH)
      CH = MOD( ABS(CH),360.0 )
      IF(SG.LE.0.0) CH=360.-CH
C
  999 CONTINUE
      END
