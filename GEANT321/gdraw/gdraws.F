*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRAWS(ISHAPE,PAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the shape number ISHAPE, of parameters PAR          *
C.    *                                                                *
C.    * SHAPE     SHAPE    SHAPE                                       *
C.    * NUMBER    TYPE     PARAMETERS                                  *
C.    * -------------------------------------------------------------- *
C.    *                                                                *
C.    *   1       BOX      DX,DY,DZ                                    *
C.    *   2       TRD1     DX1,DX2,DY,DZ                               *
C.    *   3       TRD2     DX1,DX2,DY1,DY2,DZ                          *
C.    *   4       TRAP     DZ,TX,TY,H1,BL1,TL1,TTH1,H2,BL2,TL2,TTH2    *
C.    *                                                                *
C.    *   5       TUBE     RMIN,RMAX,DZ                                *
C.    *   6       TUBS     RMIN,RMAX,DZ,PHIMIN,PHIMAX                  *
C.    *   7       CONE     DZ,RMIN1,RMAX1,RMIN2,RMAX2                  *
C.    *   8       CONS     DZ,RMIN1,RMAX1,RMIN2,RMAX2,PHIMIN,PHIMAX    *
C.    *                                                                *
C.    *   9       SPHE     RMIN,RMAX,THEMIN,THEMAX,PHIMIN,PHIMAX       *
C.    *                                                                *
C.    *  10       PARA     DX,DY,DZ,TXY,TXZ,TYZ                        *
C.    *  11       PGON     PHIMIN,DPHI,NDIV,NZ,Z(1),RMIN(1),RMAX(1),...*
C.    *  12       PCON     PHIMIN,DPHI,NZ,Z(1),RMIN(1),RMAX(1),Z(2),...*
C.    *  13       ELTU     A,B,DZ                                      *
C.    *  14       HYPE     RMIN,RMAX,DZ,PHI                            *
C.    *  NSGTRA   GTRA     DZ,TH,PHI,TWIST,Y1,XL1,XH1,TH1,Y2,XL2,XH2,..*
C.    *  NSCTUB   CTUB     RMIN,RMAX,DZ,PHIMIN,PHIMAX,LXL,LYL,LZL,LXH,.*
C.    *                                                                *
C.    *    ==>Called by : GDRAW                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *       Modification log.                                        *
C.    *        1-Jun-88 A.C.McPherson - Introduce cut tube shape.      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
*
      PARAMETER ( NLPC = 40 )
*
*            The constant NLPC defined in the parameter statement
*            is the number of line elements to form a complete
*            circle in the surface definitions for a cut tube.
*
      DIMENSION CPHIS(NLPC+1),SPHIS(NLPC+1)
*
      DIMENSION X(3,46), U(46), V(46)
      DIMENSION PAR(50),P(3,8),PP(3,8)
C.
C.    ------------------------------------------------------------------
C.
      IF (ISHAPE.NE.1) GO TO 200
C
C             BOX
C
      DX1=PAR(1)
      DY1=PAR(2)
      DX2=DX1
      DY2=DY1
      DZ=PAR(3)
      GO TO 1000
C
  200 IF (ISHAPE.NE.2) GO TO 300
C
C             TRD1
C
      DX1=PAR(1)
      DX2=PAR(2)
      DY1=PAR(3)
      DY2=DY1
      DZ=PAR(4)
      GO TO 1000
C
  300 IF (ISHAPE.NE.3) GO TO 400
C
C             TRD2
C
      DX1=PAR(1)
      DX2=PAR(2)
      DY1=PAR(3)
      DY2=PAR(4)
      DZ=PAR(5)
      GO TO 1000
C
  400 IF (ISHAPE.NE.4) GO TO 500
C
C             TRAP
C
      DZ=PAR(1)
      TX=PAR(2)
      TY=PAR(3)
      H1=PAR(4)
      BL1=PAR(5)
      TL1=PAR(6)
      TTH1=PAR(7)
      H2=PAR(8)
      BL2=PAR(9)
      TL2=PAR(10)
      TTH2=PAR(11)
      GO TO 1500
C
  500 IF (ISHAPE.NE.5) GO TO 600
C
C             TUBE
C
      RMIN1=PAR(1)
      RMAX1=PAR(2)
      RMIN2=RMIN1
      RMAX2=RMAX1
      Z2=PAR(3)
      Z1=-Z2
      GO TO 2000
C
  600 IF (ISHAPE.NE.6) GO TO 700
C
C             TUBS
C
      RMIN1=PAR(1)
      RMAX1=PAR(2)
      RMIN2=RMIN1
      RMAX2=RMAX1
      Z2=PAR(3)
      Z1=-Z2
      PHIMIN=PAR(4)
      PHIMAX=PAR(5)
      GO TO 2500
C
  700 IF (ISHAPE.NE.7) GO TO 800
C
C             CONE
C
      RMIN1=PAR(2)
      RMAX1=PAR(3)
      RMIN2=PAR(4)
      RMAX2=PAR(5)
      Z2=PAR(1)
      Z1=-Z2
      GO TO 2000
C
  800 IF (ISHAPE.NE.8) GO TO 900
C
C             CONS
C
      RMIN1=PAR(2)
      RMAX1=PAR(3)
      RMIN2=PAR(4)
      RMAX2=PAR(5)
      Z2=PAR(1)
      Z1=-Z2
      PHIMIN=PAR(6)
      PHIMAX=PAR(7)
      GO TO 2500
C
  900 IF (ISHAPE.NE.9) GO TO 910
C
C             SPHE
C
      RMIN=PAR(1)
      RMAX=PAR(2)
      PHMI=PAR(5)
      PHMA=PAR(6)
      GO TO 3000
C
  910 IF (ISHAPE.NE.10) GO TO 911
C
C             PARA
C
      DX=PAR(1)
      DY=PAR(2)
      DZ=PAR(3)
      TXY=PAR(4)
      TXZ=PAR(5)
      TYZ=PAR(6)
C
      TX=TXZ
      TY=TYZ
      H1=DY
      BL1=DX
      TL1=DX
      TTH1=TXY
      H2=DY
      BL2=DX
      TL2=DX
      TTH2=TXY
      GO TO 1500
C
  911 IF (ISHAPE.NE.11) GO TO 912
C
C             PGON
C
      PHIMIN=PAR(1)
      PHIMAX=PHIMIN+PAR(2)
      NDIV=PAR(3)
      NZ=PAR(4)
      DPHI=(PHIMAX-PHIMIN)/NDIV
C
C             Z(1)=PAR(5) - RMIN(1)=PAR(6) - RMAX(1)=PAR(7) - Z(2)=PAR(8) - ...
C
      GO TO 4000
C
  912 IF (ISHAPE.NE.12) GO TO 950
C
C             PCON
C
      PHIMIN=PAR(1)
      PHIMAX=PHIMIN+PAR(2)
      NZ=PAR(3)
C
C             Z(1)=PAR(4) - RMIN(1)=PAR(5) - RMAX(1)=PAR(6) - Z(2)=PAR(7) - ...
C
      GO TO 5000
C
  950 CONTINUE
C
      IF (ISHAPE.NE.13) GO TO 951
C
C             ELTU
C
      A=PAR(1)
      B=PAR(2)
      Z2=PAR(3)
      Z1=-Z2
      GO TO 7000
  951 CONTINUE
C
      IF (ISHAPE.NE.14) GO TO 955
C
C             HYPErboloid
C
      RMIN1 = PAR(1)
      RMAX1 = PAR(2)
      Z2=PAR(3)
      TANTHS = (TAN(PAR(4)*DEGRAD))**2
      RMIN12 = RMIN1*RMIN1
      RMAX12 = RMAX1*RMAX1
      RMIN2 = SQRT(RMIN12 + Z2*Z2*TANTHS)
      RMAX2 = SQRT(RMAX12 + Z2*Z2*TANTHS)
      Z1=-Z2
      GO TO 7410
 
  955 CONTINUE
      IF(ISHAPE.NE.28) GO TO 980
C
C             General twisted trapezoid.
C
      DO 970 IL=1,4
      I0=IL*4+11
      P(3,IL)=-PAR(1)
      P(1,IL)=PAR(I0)+PAR(I0+2)*P(3,IL)
      P(2,IL)=PAR(I0+1)+PAR(I0+3)*P(3,IL)
      P(3,IL+4)=PAR(1)
      P(1,IL+4)=PAR(I0)+PAR(I0+2)*P(3,IL+4)
      P(2,IL+4)=PAR(I0+1)+PAR(I0+3)*P(3,IL+4)
  970 CONTINUE
C
      GO TO 1600
*
  980 CONTINUE
      IF( ISHAPE .EQ. NSCTUB ) THEN
*
        DPHIS = PAR(5)-PAR(4)
        IF( DPHIS .LE. 0.0 ) DPHIS=DPHIS+TWOPI
        NL = MAX(DPHIS*NLPC/360.0,1.)
        DPHI = 360.0/NLPC
        PHIS = PAR(4)
        IF( PAR(4) .EQ. 0.0 .AND. PAR(5) .EQ. 360.0 ) THEN
          ISEG = 0
        ELSE
          ISEG = 1
          DPHI = DPHIS/NL
        ENDIF
*
        GO TO 6000
*
      ELSE
        GO TO 9999
      ENDIF
C
 1000 CONTINUE
C
C             Rectilinear shapes: BOX,TRD1,TRD2
C
      X1=0.
      Y1=0.
      X2=0.
      Y2=0.
      Z1=-DZ
      Z2=DZ
C
C             Calculate the 8 vertex for rectilinear shapes
C
      P(1,1)=X1+DX1
      P(2,1)=Y1+DY1
      P(3,1)=Z1
      P(1,2)=X1-DX1
      P(2,2)=Y1+DY1
      P(3,2)=Z1
      P(1,3)=X1-DX1
      P(2,3)=Y1-DY1
      P(3,3)=Z1
      P(1,4)=X1+DX1
      P(2,4)=Y1-DY1
      P(3,4)=Z1
      P(1,5)=X2+DX2
      P(2,5)=Y2+DY2
      P(3,5)=Z2
      P(1,6)=X2-DX2
      P(2,6)=Y2+DY2
      P(3,6)=Z2
      P(1,7)=X2-DX2
      P(2,7)=Y2-DY2
      P(3,7)=Z2
      P(1,8)=X2+DX2
      P(2,8)=Y2-DY2
      P(3,8)=Z2
C
      CALL GDRECT(P(1,1),P(1,2),P(1,3),P(1,4))
      CALL GDRECT(P(1,5),P(1,6),P(1,7),P(1,8))
      CALL GDLINE(P(1,1),P(1,5))
      CALL GDLINE(P(1,2),P(1,6))
      CALL GDLINE(P(1,3),P(1,7))
      CALL GDLINE(P(1,4),P(1,8))
C
      GO TO 9999
C
 1500 CONTINUE
C
C             TRAP,PARA
C
C             Calculate the 8 vertex
C
      P(1,1)=-DZ*TX+TTH1*H1+TL1
      P(2,1)=+H1-DZ*TY
      P(3,1)=-DZ
      P(1,2)=-DZ*TX+TTH1*H1-TL1
      P(2,2)=+H1-DZ*TY
      P(3,2)=-DZ
      P(1,3)=-DZ*TX-TTH1*H1-BL1
      P(2,3)=-H1-DZ*TY
      P(3,3)=-DZ
      P(1,4)=-DZ*TX-TTH1*H1+BL1
      P(2,4)=-H1-DZ*TY
      P(3,4)=-DZ
      P(1,5)=+DZ*TX+TTH2*H2+TL2
      P(2,5)=+H2+DZ*TY
      P(3,5)=+DZ
      P(1,6)=+DZ*TX+TTH2*H2-TL2
      P(2,6)=+H2+DZ*TY
      P(3,6)=+DZ
      P(1,7)=+DZ*TX-TTH2*H2-BL2
      P(2,7)=-H2+DZ*TY
      P(3,7)=+DZ
      P(1,8)=+DZ*TX-TTH2*H2+BL2
      P(2,8)=-H2+DZ*TY
      P(3,8)=+DZ
C
 1600 CONTINUE
C
      CALL GDRECT(P(1,1),P(1,2),P(1,3),P(1,4))
      CALL GDRECT(P(1,5),P(1,6),P(1,7),P(1,8))
      CALL GDLINE(P(1,1),P(1,5))
      CALL GDLINE(P(1,2),P(1,6))
      CALL GDLINE(P(1,3),P(1,7))
      CALL GDLINE(P(1,4),P(1,8))
C
C             Condition for plane sides are :
C
C             TTH1=TTH2
C
C             and
C
C             H2*(BL1-TL1)=H1(BL2-TL2)
C
C             In that case we should draw on each side 10 lines
C             (perpendicular to side lines) to make an easy
C             visualisation that sides are not planes
C
      GO TO 9999
C
 2000 CONTINUE
C
C             Cylindric shapes: TUBE,CONE
C
      CALL GDCIRC(RMAX1,Z1)
      CALL GDCIRC(RMIN1,Z1)
      CALL GDCIRC(RMAX2,Z2)
      CALL GDCIRC(RMIN2,Z2)
      PHIP=GPHI+90.
      PHIM=GPHI+270.
      CALL GDLCYL(RMIN1,PHIP,Z1,RMIN2,PHIP,Z2)
      CALL GDLCYL(RMIN1,PHIM,Z1,RMIN2,PHIM,Z2)
      CALL GDLCYL(RMAX1,PHIP,Z1,RMAX2,PHIP,Z2)
      CALL GDLCYL(RMAX1,PHIM,Z1,RMAX2,PHIM,Z2)
C
      GO TO 9999
C
 2500 CONTINUE
C
C             Segmented cylindric shapes: TUBS,CONS
C
      CALL GDARC(RMAX1,Z1,PHIMIN,PHIMAX)
      CALL GDARC(RMIN1,Z1,PHIMIN,PHIMAX)
      CALL GDARC(RMAX2,Z2,PHIMIN,PHIMAX)
      CALL GDARC(RMIN2,Z2,PHIMIN,PHIMAX)
      PHIP=AMOD((GPHI+90.),360.)
      PHIM=AMOD((GPHI+270.),360.)
      IF (PHIP.LE.PHIMIN.OR.PHIP.GE.PHIMAX) GO TO 2510
      CALL GDLCYL(RMIN1,PHIP,Z1,RMIN2,PHIP,Z2)
      CALL GDLCYL(RMAX1,PHIP,Z1,RMAX2,PHIP,Z2)
 2510 IF (PHIM.LE.PHIMIN.OR.PHIM.GE.PHIMAX) GO TO 2520
      CALL GDLCYL(RMIN1,PHIM,Z1,RMIN2,PHIM,Z2)
      CALL GDLCYL(RMAX1,PHIM,Z1,RMAX2,PHIM,Z2)
 2520 CALL GDLCYL(RMAX1,PHIMIN,Z1,RMAX2,PHIMIN,Z2)
      CALL GDLCYL(RMAX1,PHIMAX,Z1,RMAX2,PHIMAX,Z2)
      CALL GDLCYL(RMIN1,PHIMIN,Z1,RMIN2,PHIMIN,Z2)
      CALL GDLCYL(RMIN1,PHIMAX,Z1,RMIN2,PHIMAX,Z2)
      CALL GDLCYL(RMAX1,PHIMIN,Z1,RMIN1,PHIMIN,Z1)
      CALL GDLCYL(RMAX2,PHIMIN,Z2,RMIN2,PHIMIN,Z2)
      CALL GDLCYL(RMAX1,PHIMAX,Z1,RMIN1,PHIMAX,Z1)
      CALL GDLCYL(RMAX2,PHIMAX,Z2,RMIN2,PHIMAX,Z2)
C
      GO TO 9999
C
 3000 CONTINUE
C
C             SPHE
C
      CALL GDARC(RMAX,0.,PHMI,PHMA)
      CALL GDARC(RMIN,0.,PHMI,PHMA)
      DP = PHMA-PHMI
      IF(DP.LE.0.) DP=DP+360.
      NSTEP = MAX(DP/15.,1.)
      DDP = DP/NSTEP
      PPH = PHMI-DDP
      DO 3005 I=1,NSTEP+1
      PPH = PPH+DDP
      RPPH = PPH*DEGRAD
      COSPH = COS(RPPH)
      SINPH = SIN(RPPH)
      DO 3004 J=1,46
      THET = (J-1)*PI/45.
      X(1,J) = RMAX*SIN(THET)*COSPH
      X(2,J) = RMAX*SIN(THET)*SINPH
      X(3,J) = RMAX*COS(THET)
 3004 CONTINUE
      CALL GDFR3D(X,46,U,V)
      CALL GDRAWV(U,V,46)
 3005 CONTINUE
      IF(RMIN.GE.0.) THEN
      PPH = PHMI-DDP
      DO 3007 I=1,NSTEP+1
      PPH = PPH+DDP
      RPPH = PPH*DEGRAD
      COSPH = COS(RPPH)
      SINPH = SIN(RPPH)
      DO 3006 J=1,46
      THET = (J-1)*PI/45.
      X(1,J) = RMIN*SIN(THET)*COSPH
      X(2,J) = RMIN*SIN(THET)*SINPH
      X(3,J) = RMIN*COS(THET)
 3006 CONTINUE
      CALL GDFR3D(X,46,U,V)
      CALL GDRAWV(U,V,46)
 3007 CONTINUE
      ENDIF
      DO 3010 I=1,3
      DO 3010 J=1,6
 3010 P(I,J)=0.
      IF(DP.GE.360.) THEN
      P(3,1)=-RMAX
      P(3,2)=RMAX
      P(1,3)=-RMAX
      P(1,4)=RMAX
      P(2,5)=RMAX
      P(2,6)=-RMAX
      CALL GDLINE(P(1,1),P(1,2))
      CALL GDLINE(P(1,3),P(1,4))
      CALL GDLINE(P(1,5),P(1,6))
      ELSE
      P(1,1) = RMIN*COS(PHMI*DEGRAD)
      P(2,1) = RMIN*SIN(PHMI*DEGRAD)
      P(1,2) = RMAX*COS(PHMI*DEGRAD)
      P(2,2) = RMAX*SIN(PHMI*DEGRAD)
      CALL GDLINE(P(1,1),P(1,2))
      P(1,1) = RMIN*COS(PHMA*DEGRAD)
      P(2,1) = RMIN*SIN(PHMA*DEGRAD)
      P(1,2) = RMAX*COS(PHMA*DEGRAD)
      P(2,2) = RMAX*SIN(PHMA*DEGRAD)
      CALL GDLINE(P(1,1),P(1,2))
      P(3,3) = -RMAX
      P(3,4) = -RMIN
      CALL GDLINE(P(1,3),P(1,4))
      P(3,3) =  RMAX
      P(3,4) =  RMIN
      CALL GDLINE(P(1,3),P(1,4))
      ENDIF
      GO TO 9999
C
 4000 CONTINUE
C
C             PGON
C
      FACT=1./COS(DEGRAD*DPHI/2.)
      DO 4002 IZ=1,NZ
        PAR(6+(IZ-1)*3)=PAR(6+(IZ-1)*3)*FACT
        PAR(7+(IZ-1)*3)=PAR(7+(IZ-1)*3)*FACT
 4002 CONTINUE
C
      DO 4050 IZ=1,NZ
C
        ZI=PAR(5+(IZ-1)*3)
        R0=PAR(6+(IZ-1)*3)
        R1=PAR(7+(IZ-1)*3)
C
        IF (IZ.EQ.1.OR.IZ.EQ.NZ) GO TO 4003
        R0PRE=PAR(6+(IZ-2)*3)
        R0POST=PAR(6+IZ*3)
        IF (R0.EQ.R0PRE)GO TO 4006
        IF (R0.EQ.R0POST)GO TO 4006
 4003     CONTINUE
          DO 4005 IDIV=1,NDIV
            PHI0=PHIMIN+(IDIV-1)*DPHI
            PHI1=PHI0+DPHI
            CALL GDLCYL(R0,PHI0,ZI,R0,PHI1,ZI)
 4005     CONTINUE
C
 4006   IF (IZ.EQ.1.OR.IZ.EQ.NZ) GO TO 4008
        R1PRE=PAR(7+(IZ-2)*3)
        R1POST=PAR(7+IZ*3)
        IF (R1.EQ.R1PRE )GO TO 4020
        IF (R1.EQ.R1POST)GO TO 4020
 4008     CONTINUE
          DO 4010 IDIV=1,NDIV
            PHI0=PHIMIN+(IDIV-1)*DPHI
            PHI1=PHI0+DPHI
            CALL GDLCYL(R1,PHI0,ZI,R1,PHI1,ZI)
 4010     CONTINUE
C
 4020   IF ((IZ.EQ.1.OR.IZ.EQ.NZ).AND.(PHIMAX-PHIMIN.NE.360.)) THEN
          CALL GDLCYL(R0,PHIMIN,ZI,R1,PHIMIN,ZI)
          CALL GDLCYL(R0,PHIMAX,ZI,R1,PHIMAX,ZI)
        ENDIF
C
        IF (IZ.EQ.1) GO TO 4050
C
        ZI0=PAR(5+(IZ-2)*3)
        R00=PAR(6+(IZ-2)*3)
        R10=PAR(7+(IZ-2)*3)
        DO 4030 IDIV=1,NDIV
          PH=PHIMIN+(IDIV-1)*DPHI
          CALL GDLCYL(R00,PH,ZI0,R0,PH,ZI)
          CALL GDLCYL(R10,PH,ZI0,R1,PH,ZI)
 4030   CONTINUE
        CALL GDLCYL(R00,PHIMAX,ZI0,R0,PHIMAX,ZI)
        CALL GDLCYL(R10,PHIMAX,ZI0,R1,PHIMAX,ZI)
C
 4050 CONTINUE
C
      GO TO 9999
C
 5000 CONTINUE
C
C             PCON
C
      DO 5555 IZ=1,NZ
        ZI=PAR(4+(IZ-1)*3)
        R0=PAR(5+(IZ-1)*3)
        R1=PAR(6+(IZ-1)*3)
        IF (IZ.EQ.1.OR.IZ.EQ.NZ) GO TO 5010
        R1PRE=PAR(6+(IZ-2)*3)
        R1POST=PAR(6+IZ*3)
        IF (R1.LE.R1PRE.OR.R1.LE.R1POST) GO TO 5015
 5010   CONTINUE
        CALL GDARC(R0,ZI,PHIMIN,PHIMAX)
        CALL GDARC(R1,ZI,PHIMIN,PHIMAX)
 5015   CONTINUE
        IF ((PHIMAX-PHIMIN).EQ.360.) GO TO 5020
        CALL GDLCYL(R0,PHIMIN,ZI,R1,PHIMIN,ZI)
        CALL GDLCYL(R0,PHIMAX,ZI,R1,PHIMAX,ZI)
 5020   CONTINUE
        IF (IZ.EQ.1) GO TO 5555
        ZI0=PAR(4+(IZ-2)*3)
        R00=PAR(5+(IZ-2)*3)
        R10=PAR(6+(IZ-2)*3)
        IF ((PHIMAX-PHIMIN).EQ.360.) GO TO 5030
        CALL GDLCYL(R00,PHIMIN,ZI0,R0,PHIMIN,ZI)
        CALL GDLCYL(R10,PHIMIN,ZI0,R1,PHIMIN,ZI)
        CALL GDLCYL(R00,PHIMAX,ZI0,R0,PHIMAX,ZI)
        CALL GDLCYL(R10,PHIMAX,ZI0,R1,PHIMAX,ZI)
C
 5030   CONTINUE
        PHIP=AMOD((GPHI+90.),360.)
        PHIM=AMOD((GPHI+270.),360.)
        IF (PHIP.LT.PHIMIN.OR.PHIP.GT.PHIMAX) GO TO 5510
        CALL GDLCYL(R00,PHIP,ZI0,R0,PHIP,ZI)
        CALL GDLCYL(R10,PHIP,ZI0,R1,PHIP,ZI)
 5510   IF (PHIM.LT.PHIMIN.OR.PHIM.GT.PHIMAX) GO TO 5555
        CALL GDLCYL(R00,PHIM,ZI0,R0,PHIM,ZI)
        CALL GDLCYL(R10,PHIM,ZI0,R1,PHIM,ZI)
 5555 CONTINUE
C
      GO TO 9999
*
 6000 CONTINUE
*
*          Cut tube shape.
*
      CPHIS(1) = COS( PHIS*DEGRAD )
      SPHIS(1) = SIN( PHIS*DEGRAD )
      DO 6010 I = 1, NL
        PHIS = PHIS+DPHI
        CPHIS(I+1) = COS( PHIS*DEGRAD )
        SPHIS(I+1) = SIN( PHIS*DEGRAD )
 6010 CONTINUE
      P( 1, 1) = PAR(2)*CPHIS(1)
      P( 2, 1) = PAR(2)*SPHIS(1)
      P( 3, 1) = -( PAR( 6)*P( 1, 1) +PAR( 7)*P( 2, 1) )/PAR( 8)
     + -PAR( 3)
      P( 1, 4) = PAR(1)*CPHIS(1)
      P( 2, 4) = PAR(1)*SPHIS(1)
      P( 3, 4) = -( PAR( 6)*P( 1, 4) +PAR( 7)*P( 2, 4) )/PAR( 8)
     + -PAR( 3)
      P( 1, 5) = PAR(2)*CPHIS(1)
      P( 2, 5) = PAR(2)*SPHIS(1)
      P( 3, 5) = -( PAR( 9)*P( 1, 5) +PAR( 10)*P( 2, 5) )/PAR( 11)
     + +PAR( 3)
      P( 1, 8) = PAR(1)*CPHIS(1)
      P( 2, 8) = PAR(1)*SPHIS(1)
      P( 3, 8) = -( PAR( 9)*P( 1, 8) +PAR( 10)*P( 2, 8) )/PAR( 11)
     + +PAR( 3)
*
      IF( ISEG .EQ. 1 ) THEN
        CALL GDRECT( P( 1, 1), P( 1, 4), P( 1, 8), P( 1, 5) )
      ENDIF
*
      DO 6020 I = 1, NL
        P( 1, 2) = PAR(2)*CPHIS(I+1)
        P( 2, 2) = PAR(2)*SPHIS(I+1)
        P( 3, 2) = -( PAR( 6)*P( 1, 2) +PAR( 7)*P( 2, 2) )/PAR( 8)
     +   -PAR( 3)
        CALL GDLINE( P( 1, 1),  P( 1, 2) )
        P( 1, 1) = P( 1, 2)
        P( 2, 1) = P( 2, 2)
        P( 3, 1) = P( 3, 2)
 6020 CONTINUE
*
      DO 6030 I = 1, NL
        P( 1, 3) = PAR(1)*CPHIS(I+1)
        P( 2, 3) = PAR(1)*SPHIS(I+1)
        P( 3, 3) = -( PAR( 6)*P( 1, 3) +PAR( 7)*P( 2, 3) )/PAR( 8)
     +   -PAR( 3)
        CALL GDLINE( P( 1, 4),  P( 1, 3) )
        P( 1, 4) = P( 1, 3)
        P( 2, 4) = P( 2, 3)
        P( 3, 4) = P( 3, 3)
 6030 CONTINUE
*
      DO 6040 I = 1, NL
        P( 1, 6) = PAR(2)*CPHIS(I+1)
        P( 2, 6) = PAR(2)*SPHIS(I+1)
        P( 3, 6) = -( PAR( 9)*P( 1, 6) +PAR( 10)*P( 2, 6) )/PAR( 11)
     +   +PAR( 3)
        CALL GDLINE( P( 1, 5),  P( 1, 6) )
        P( 1, 5) = P( 1, 6)
        P( 2, 5) = P( 2, 6)
        P( 3, 5) = P( 3, 6)
 6040 CONTINUE
*
      DO 6050 I = 1, NL
        P( 1, 7) = PAR(1)*CPHIS(I+1)
        P( 2, 7) = PAR(1)*SPHIS(I+1)
        P( 3, 7) = -( PAR( 9)*P( 1, 7) +PAR( 10)*P( 2, 7) )/PAR( 11)
     +   +PAR( 3)
        CALL GDLINE( P( 1, 8),  P( 1, 7) )
        P( 1, 8) = P( 1, 7)
        P( 2, 8) = P( 2, 7)
        P( 3, 8) = P( 3, 7)
 6050 CONTINUE
*
      IF( ISEG .EQ. 1 ) THEN
        CALL GDRECT( P( 1, 1), P( 1, 4), P( 1, 8), P( 1, 5) )
      ENDIF
*
      PHIP = AMOD( GPHI+90.0, 360.0 )
      PHIM = AMOD( GPHI+270.0, 360.0 )
      DPHIP = PHIP-PAR(4)
      DPHIM = PHIM-PAR(4)
      IF( DPHIP .LT. 0.0 ) DPHIP = DPHIP+TWOPI
      IF( DPHIM .LT. 0.0 ) DPHIM = DPHIM+TWOPI
*
      IF( DPHIP .LE. DPHIS ) THEN
        CP = COS( PHIP*DEGRAD )
        SP = SIN( PHIP*DEGRAD )
        P( 1, 1) = PAR(2)*CP
        P( 2, 1) = PAR(2)*SP
        P( 3, 1) = -( PAR( 6)*P( 1, 1) +PAR( 7)*P( 2, 1) )/PAR( 8)
     +   -PAR( 3)
        P( 1, 4) = PAR(1)*CP
        P( 2, 4) = PAR(1)*SP
        P( 3, 4) = -( PAR( 6)*P( 1, 4) +PAR( 7)*P( 2, 4) )/PAR( 8)
     +   -PAR( 3)
        P( 1, 5) = PAR(2)*CP
        P( 2, 5) = PAR(2)*SP
        P( 3, 5) = -( PAR( 9)*P( 1, 5) +PAR( 10)*P( 2, 5) )/PAR( 11)
     +   +PAR( 3)
        P( 1, 8) = PAR(1)*CP
        P( 2, 8) = PAR(1)*SP
        P( 3, 8) = -( PAR( 9)*P( 1, 8) +PAR( 10)*P( 2, 8) )/PAR( 11)
     +   +PAR( 3)
        CALL GDLINE( P( 1, 1), P( 1, 5) )
        CALL GDLINE( P( 1, 4), P( 1, 8) )
*
      ENDIF
*
      IF( DPHIM .LE. DPHIS ) THEN
        CP = COS( PHIM*DEGRAD )
        SP = SIN( PHIM*DEGRAD )
        P( 1, 1) = PAR(2)*CP
        P( 2, 1) = PAR(2)*SP
        P( 3, 1) = -( PAR( 6)*P( 1, 1) +PAR( 7)*P( 2, 1) )/PAR( 8)
     +   -PAR( 3)
        P( 1, 4) = PAR(1)*CP
        P( 2, 4) = PAR(1)*SP
        P( 3, 4) = -( PAR( 6)*P( 1, 4) +PAR( 7)*P( 2, 4) )/PAR( 8)
     +   -PAR( 3)
        P( 1, 5) = PAR(2)*CP
        P( 2, 5) = PAR(2)*SP
        P( 3, 5) = -( PAR( 9)*P( 1, 5) +PAR( 10)*P( 2, 5) )/PAR( 11)
     +   +PAR( 3)
        P( 1, 8) = PAR(1)*CP
        P( 2, 8) = PAR(1)*SP
        P( 3, 8) = -( PAR( 9)*P( 1, 8) +PAR( 10)*P( 2, 8) )/PAR( 11)
     +   +PAR( 3)
        CALL GDLINE( P( 1, 1), P( 1, 5) )
        CALL GDLINE( P( 1, 4), P( 1, 8) )
*
      ENDIF
      GO TO 9999
C
 7000 CONTINUE
C
C             ELTU
C
      CALL GDELTU(A,B,Z1)
      CALL GDELTU(A,B,Z2)
      P(1,1)=A
      P(2,1)=0.
      P(3,1)=Z1
      P(1,2)=A
      P(2,2)=0.
      P(3,2)=Z2
      P(1,3)=-A
      P(2,3)=0.
      P(3,3)=Z1
      P(1,4)=-A
      P(2,4)=0.
      P(3,4)=Z2
      CALL GDLINE(P(1,1),P(1,2))
      CALL GDLINE(P(1,3),P(1,4))
      GO TO 9999
C                               draw HYPErboloid
 7410 CONTINUE
      PP(2,1) = RMAX2
      PP(2,3) = RMIN2
      P(1,1) = 0.
      P(1,2) = 0.
      P(1,3) = 0.
      P(1,4) = 0.
      PP(3,1) = Z2
      PP(3,3) = Z2
      CALL GDCIRC(RMAX2,Z1)
      CALL GDCIRC(RMIN2,Z1)
      CALL GDCIRC(RMAX2,Z2)
      CALL GDCIRC(RMIN2,Z2)
      NZSTEP = 20
      DELZ = Z2 / NZSTEP
      DO 7440 IZ = 1, NZSTEP
         ZZ = Z2 - IZ*DELZ
         PP(3,2) = ZZ
         PP(3,4) = ZZ
         ZZZ = ZZ*ZZ*TANTHS
         PP(2,2) = SQRT(RMAX12 + ZZZ)
         PP(2,4) = SQRT(RMIN12 + ZZZ)
         DO 7430 ISY = -1, +1, 2
         DO 7430 ISZ = -1, +1, 2
            DO 7420 J = 1, 4
               P(2,J) = ISY * PP(2,J)
               P(3,J) = ISZ * PP(3,J)
 7420       CONTINUE
            CALL GDLINE(P(1,1),P(1,2))
            CALL GDLINE(P(1,3),P(1,4))
 7430    CONTINUE
         PP(2,1) = PP(2,2)
         PP(2,3) = PP(2,4)
         PP(3,1) = PP(3,2)
         PP(3,3) = PP(3,4)
 7440 CONTINUE
C
 9999 END
