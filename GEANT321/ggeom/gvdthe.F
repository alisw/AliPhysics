*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GVDTHE(ISH,IROT,DX,PARS,CL,CH,IERR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO COMPUTE THE THETA LIMITS FOR VOLUME OF SHAPE     *
C.    *    ISH ROTATED BY MATRIX IROT AND DISPLACED BY VECTOR DX.      *
C.    *    THE VOLUME HAS NPAR PARAMETERS IN THE ARRAY PARS. THE LOWER *
C.    *    LIMIT IS RETURNED IN CL THE HIGHER IN CH. IF THE            *
C.    *    CALCULATION CANNOT BE MADE IERR IS SET TO 1 OTHERWISE IT    *
C.    *    IS SET TO 0.                                                *
C.    *                                                                *
C.    *    ==>Called by : GVDLIM                                       *
C.    *         Author  S.Giani  *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
      DIMENSION DX(3),PARS(50),X(3),XT(3),X1(3),X2(3),XT1(3),XT2(3)
C.
C.          ----------------------------------------------
C.
      IERR=1
C
      DXS=DX(1)*DX(1)+DX(2)*DX(2)
      IF(DXS.GT.0.0) DXS=SQRT(DXS)
C
      IF(ISH.GT.4.AND.ISH.NE.10.AND.ISH.NE.28) GO TO 130
C
C            CUBES, TRAPEZOIDS AND PARALLELEPIPEDS.
C
      IERR=0
      CL=0.0
      CH=180.0
C
      IF(DXS.LE.0.0) GO TO 999
C
      TH=90.
      IF(ABS(DX(3)).LT.1.0E-06)GO TO 5
      TH=ATAN(DXS/DX(3))*RADDEG
      IF(TH.LT.0.0) TH=TH+180.0
   5  TL=TH
C
      DO 50 IP=1,8
C
C           THIS IS A LOOP OVER THE 8 CORNERS.
C           FIRST FIND THE LOCAL COORDINATES.
C
      IF(ISH.EQ.28) THEN
         goto 999
C            General twisted trapezoid.
C
*         IL=(IP+1)/2
*         I0=IL*4+11
*         IS=(IP-IL*2)*2+1
*         X(3)=PARS(1)*IS
*         X(1)=PARS(I0)+PARS(I0+2)*X(3)
*         X(2)=PARS(I0+1)+PARS(I0+3)*X(3)
*         GO TO 20
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
      Z=DX(3)+XT(3)
      R=(DX(1)+XT(1))**2+(DX(2)+XT(2))**2
      RT=R+(DX(3)+XT(3))**2
      IF(RT.LT.1.0E-10) GO TO 999
C
      IF(R.GT.0.0) GO TO 30
      IF(Z.GT.0.0) T=0.0
      IF(Z.LT.0.0) T=180.0
      GO TO 40
C
   30 CONTINUE
      T=90.0
      IF(ABS(Z).LT.1.0E-6) GO TO 40
C
      R=SQRT(R)
      T=ATAN(R/Z)*RADDEG
      IF(T.EQ.0.0.AND.Z.LT.0.0) T=180.0
      IF(T.LT.0.0) T=T+180.0
C
   40 CONTINUE
      IF(T.GT.TH) TH=T
      IF(T.LT.TL) TL=T
C
   50 CONTINUE
C
C           THETA LIMITS SET FROM THE POINTS NOW DO THE EDGES.
C
      DO 120 IL=1,12
C
C           FIND THE END POINT NUMBERS FOR EACH EDGE.
C
      IF(IL.GT.4) GO TO 60
      IPP1=1
      IF(IL.GT.2) IPP1=4
      IPP2=2
      IF(MOD(IL,2).EQ.1) IPP2=3
C
      GO TO 80
   60 CONTINUE
      IF(IL.LT.9) GO TO 70
C
      IPP1=5
      IF(IL.GT.10) IPP1=8
      IPP2=6
      IF(MOD(IL,2).EQ.1) IPP2=7
C
      GO TO 80
   70 CONTINUE
C
      IPP1=IL-4
      IPP2=IL
C
   80 CONTINUE
C
C           NOW GET THE POINTS AND ROTATE THEM.
C
      IF(ISH.EQ.28) THEN
C
C            General twisted trapezoid.
C
         ILL=IPP1
         IF(IPP1.EQ.3) ILL=4
         IF(IPP1.EQ.4) ILL=3
         I0=ILL*4+11
         X1(3)=PARS(1)
         IF(IPP1.LT.5) X1(3)=-X1(3)
         X1(1)=PARS(I0)+PARS(I0+2)*X1(3)
         X1(2)=PARS(I0+1)+PARS(I0+3)*X1(3)
         ILL=IPP2
         IF(IPP2.EQ.3) ILL=4
         IF(IPP2.EQ.4) ILL=3
         I0=ILL*4+11
         X2(3)=PARS(1)
         IF(IPP2.LT.5) X2(3)=-X2(3)
         X2(1)=PARS(I0)+PARS(I0+2)*X2(3)
         X2(2)=PARS(I0+1)+PARS(I0+3)*X2(3)
C
         GO TO 100
C
      ENDIF
C
      IP3=ISH+2
      IF(ISH.EQ.10) IP3=3
      IF(ISH.EQ.4) IP3=1
      X1(3)=PARS(IP3)
      IF(IPP1.LE.4) X1(3)=-X1(3)
      IP2=3
      IF(ISH.GT.2.AND.X1(3).GT.0.0) IP2=4
      IF(ISH.EQ.1.OR.ISH.EQ.10) IP2=2
      IF(ISH.EQ.4) IP2=4
      IF(ISH.EQ.4.AND.X1(3).GT.0.0) IP2=8
      X1(2)=PARS(IP2)
      IF(MOD(IPP1+3,4).LT.2) X1(2)=-X1(2)
      IP1=1
      IF(ISH.NE.1.AND.ISH.NE.10.AND.X1(3).GT.0.0) IP1=2
      IF(ISH.EQ.4) IP1=5
      IF(ISH.EQ.4.AND.X1(3).GT.0.0) IP1=IP1+4
      IF(ISH.EQ.4.AND.X1(2).GT.0.0) IP1=IP1+1
      X1(1)=PARS(IP1)
      IF(MOD(IPP1,2).EQ.1) X1(1)=-X1(1)
C
      IP3=ISH+2
      IF(ISH.EQ.10) IP3=3
      IF(ISH.EQ.4) IP3=1
      X2(3)=PARS(IP3)
      IF(IPP2.LE.4) X2(3)=-X2(3)
      IP2=3
      IF(ISH.GT.2.AND.X2(3).GT.0.0) IP2=4
      IF(ISH.EQ.1.OR.ISH.EQ.10) IP2=2
      IF(ISH.EQ.4) IP2=4
      IF(ISH.EQ.4.AND.X2(3).GT.0.0) IP2=8
      X2(2)=PARS(IP2)
      IF(MOD(IPP2+3,4).LT.2) X2(2)=-X2(2)
      IP1=1
      IF(ISH.NE.1.AND.ISH.NE.10.AND.X2(3).GT.0.0) IP1=2
      IF(ISH.EQ.4) IP1=5
      IF(ISH.EQ.4.AND.X2(3).GT.0.0) IP1=IP1+4
      IF(ISH.EQ.4.AND.X2(2).GT.0.0) IP1=IP1+1
      X2(1)=PARS(IP1)
      IF(MOD(IPP2,2).EQ.1) X2(1)=-X2(1)
C
      IF(ISH.NE.10) GO TO 90
      X1(1)=X1(1)+X1(2)*PARS(4)+X1(3)*PARS(5)
      X1(2)=X1(2)+X1(3)*PARS(6)
      X2(1)=X2(1)+X2(2)*PARS(4)+X2(3)*PARS(5)
      X2(2)=X2(2)+X2(3)*PARS(6)
   90 CONTINUE
C
      IF(ISH.NE.4) GO TO 100
      IP4=7
      IF(X1(3).GT.0.0) IP4=11
      X1(1)=X1(1)+X1(2)*PARS(IP4)+X1(3)*PARS(2)
      X1(2)=X1(2)+X1(3)*PARS(3)
      IP4=7
      IF(X2(3).GT.0.0) IP4=11
      X2(1)=X2(1)+X2(2)*PARS(IP4)+X2(3)*PARS(2)
      X2(2)=X2(2)+X2(3)*PARS(3)
  100 CONTINUE
C
C          ROTATE.
C
      JROT=LQ(JROTM-IROT)
      XT1(1)=X1(1)
      XT1(2)=X1(2)
      XT1(3)=X1(3)
      IF(IROT.NE.0) CALL GINROT(X1,Q(JROT+1),XT1)
      XT2(1)=X2(1)
      XT2(2)=X2(2)
      XT2(3)=X2(3)
      IF(IROT.NE.0) CALL GINROT(X2,Q(JROT+1),XT2)
C
C           NOW WE HAVE THE END POINTS IN THE MASTER SYSTEM.
C           FIND THE TANGENT POINT TO THE SET OF CONES OF CONSTANT
C           THETA.
C
      DS=(XT2(1)-XT1(1))**2+(XT2(2)-XT1(2))**2+(XT2(3)-XT1(3))**2
      IF(DS.LE.0.0) GO TO 120
C
      DS=SQRT(DS)
C
      X0=(XT2(1)+XT1(1))/2.0+DX(1)
      Y0=(XT2(2)+XT1(2))/2.0+DX(2)
      Z0=(XT2(3)+XT1(3))/2.0+DX(3)
      ALX=(XT2(1)-XT1(1))/DS
      ALY=(XT2(2)-XT1(2))/DS
      ALZ=(XT2(3)-XT1(3))/DS
C
      T=90.0
      IF(Z0.EQ.0.0.AND.ALZ.EQ.0.0) GO TO 110
C
      IF(ALX.EQ.0.0.AND.ALY.EQ.0.0) GO TO 120
C             THIS CHECKS WHETHER THE LINE IS PARALLEL TO THE
C             Z AXIS IN WHICH CASE THERE IS NO TANGENT AND
C             THE END POINTS DETERMINE THE THETA RANGE.
C
      SNUM=(X0*Z0*ALX+Y0*Z0*ALY+X0*X0*ALZ+Y0*Y0*ALZ)
      SDEN=(Z0*ALX*ALX-X0*ALX*ALZ+Z0*ALY*ALY-Y0*ALY*ALZ)
C
      IF(ABS(SNUM).GT.0.5*DS*ABS(SDEN)) GO TO 120
C
C           TANGENT EXIST BETWEEN THE TWO ENDS.
C
      S = SNUM/SDEN
      X0=X0+S*ALX
      Y0=Y0+S*ALY
      Z0=Z0+S*ALZ
C
      R=X0*X0+Y0*Y0
      RT=R+Z0*Z0
C
      IF(RT.LT.1.0E-10) GO TO 999
C
      IF(R.GT.0.0) R=SQRT(R)
      T=90.0
      IF(ABS(Z0).LT.1.0E-06) GO TO 110
C
      T=ATAN(R/Z0)*RADDEG
      IF(T.EQ.0.0.AND.Z0.LT.0.0) T=180.0
      IF(T.LT.0.0) T=T+180.0
C
  110 CONTINUE
C
      IF(T.LT.TL) TL=T
      IF(T.GT.TH) TH=T
C
C        CHECK FOR THE POSSIBILITY OF STRADDLING T=0.0 &/OR 180.0.
C
      C=X0*DX(1)+Y0*DX(2)
      IF(C.GT.0.0) GO TO 120
C
C          CHECK IF SAME SIGN OF Z.
C
      IF(Z0*DX(3).LT.0.0) GO TO 999
C
      T=0.0
      IF(Z0.LT.0.0) T=180.0
      IF(T.LT.TL) TL=T
      IF(T.GT.TH) TH=T
C
  120 CONTINUE
C
C            DONE SET THE LIMITS.
C
      CL=TL
      CH=TH
C
      GO TO 999
C
  130 CONTINUE
C
C          TUBES, CONES ETC.
C
      IF(IROT.NE.0.OR.DXS.GT.1.0E-05) GO TO 180
C
C          UNROTATED AND CENTERED.
C
      IF(ISH.GT.8.AND.ISH.NE.NSCTUB.AND.ISH.NE.13
     +           .AND.ISH.NE.14)GO TO 170
C
C               TUBES AND CONES.
C
      IERR=0
      DZ=PARS(3)
      RMN=PARS(1)
      RMX=PARS(2)
      IF(ISH.EQ.13) THEN
**
**         approxime to a cylinder whit radius
**         equal to the ellipse major axis
**
        RMN=0.0
        IF(PARS(1).GT.RMX) RMX=PARS(1)
        GOTO 140
      ENDIF
      IF(ISH.EQ.14) THEN
C   not really sure of the function of these... keep RMN=PARS(1)
        RMX = SQRT(PARS(2)**2+(PARS(3)*TAN(PARS(4)*DEGRAD))**2)
        GO TO 140
      ENDIF
      IF(ISH.EQ.NSCTUB) THEN
        S1 = (1.0-PARS(8))*(1.0+PARS(8))
        IF( S1 .GT. 0.0) S1 = SQRT(S1)
        S2 = (1.0-PARS(11))*(1.0+PARS(11))
        IF( S2 .GT. 0.0) S2 = SQRT(S2)
        IF( S2 .GT. S1 ) S1 = S2
        DZ = DZ+RMX*S1
      ENDIF
      IF(ISH.LE.6.OR.ISH.EQ.NSCTUB) GO TO 140
      DZ=PARS(1)
      RMN=PARS(2)
      RMX=PARS(3)
      IF(PARS(4).LT.RMN) RMN=PARS(4)
      IF(PARS(5).GT.RMX) RMX=PARS(5)
  140 CONTINUE
C
      IF(DZ.GT.ABS(DX(3))) GO TO 160
C
C          ALL FORWARD OR ALL BACK.
C
      DSH=DX(3)-DZ
      DLN=DX(3)+DZ
      IF(DX(3).GT.0.0) GO TO 150
      DSS=DSH
      DSH=DLN
      DLN=DSS
  150 CONTINUE
C
      CL=90.0
      CH=90.0
      IF(DLN.NE.0.0) CL=ATAN(RMN/DLN)*RADDEG
      IF(DSH.NE.0.0) CH=ATAN(RMX/DSH)*RADDEG
      IF(DX(3).GT.0.0) GO TO 999
      CS=CL
      CL=CH
      CH=CS
      IF(CH.EQ.0.0) CH=180.0
      IF(CL.LT.0.0) CL=CL+180.0
      IF(CH.LT.0.0) CH=CH+180.0
C
      GO TO 999
C
  160 CONTINUE
C
      CL=90.0
      CH=90.0
      IF(DZ+DX(3).NE.0.0) CL=ATAN(RMN/(DZ+DX(3)))*RADDEG
      IF(-DZ+DX(3).NE.0.0) CH=ATAN(RMN/(-DZ+DX(3)))*RADDEG
      IF(CH.LE.0.0) CH=CH+180.0
C
      GO TO 999
C
  170 CONTINUE
      IF(ISH.GT.9) GO TO 999
C
C           SPHERE.
C
      IERR=0
      CL=PARS(3)
      CH=PARS(4)
C
      GO TO 999
  180 CONTINUE
C
      IF(ISH.EQ.11.OR.ISH.EQ.12) GOTO 999
**
      RM=PARS(2)
      IF(ISH.EQ.9) GO TO 200
C
      DZ=PARS(3)
      IF(ISH.EQ.NSCTUB) THEN
        S1 = (1.0-PARS(8))*(1.0+PARS(8))
        IF( S1 .GT. 0.0) S1 = SQRT(S1)
        S2 = (1.0-PARS(11))*(1.0+PARS(11))
        IF( S2 .GT. 0.0) S2 = SQRT(S2)
        IF( S2 .GT. S1 ) S1 = S2
        DZ = DZ+RM*S1
      ENDIF
**
      IF(ISH.EQ.13) THEN
        IF(PARS(1).GT.RM) RM=PARS(1)
        GOTO 190
      ENDIF
**
      IF(ISH.LE.6.OR.ISH.EQ.NSCTUB) GO TO 190
C
      DZ=PARS(1)
      RM=PARS(3)
      IF(PARS(5).GT.RM) RM=PARS(5)
C
  190 CONTINUE
C
      RM=SQRT(RM**2+DZ**2)
C
  200 CONTINUE
C
      CL=0.0
      CH=180.0
      IERR=0
      RC=DXS**2+DX(3)**2
      IF(RC.GT.0.0) RC=SQRT(RC)
      IF(RM.GE.RC) GO TO 999
C
      TC=90.0
      IF(ABS(DX(3)).GT.0.0) TC=ATAN(DXS/DX(3))*RADDEG
      IF(TC.LT.0.0) TC=TC+180.0
C
      DT=ASIN(RM/RC)*RADDEG
      CL=TC-DT
      CH=TC+DT
      IF(CL.LT.0.0) CL=0.0
      IF(CH.GT.180.0) CH=180.0
C
  999 CONTINUE
      END
