*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRWSC(ISHAPE,PAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the shape number ISHAPE, of parameters PAR,         *
C.    *       in cut-mode.                                             *
C.    *                                                                *
C.    *    ==>Called by : GDRAW                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *       Modification log.                                        *
C.    *       31-May-88 A.C.McPherson - Introduce cut tube shape.      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
      PARAMETER ( NLPC = 40 )
*
*            The constant NLPC defined in the parameter statement
*            is the number of line elements to form a complete
*            circle in the surface definitions for a cut tube.
*
      DIMENSION PAR(50),P(3,8)
      DIMENSION XMAX1(2,42),XMIN1(2,42),XMAX2(2,42),XMIN2(2,42)
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
      GO TO 7350
 
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
        NL = DPHIS*NLPC/360.0
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
 1600 CONTINUE
C
C             Store all the surfaces (back,front,top,bottom,right,left)
C
      CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
      CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
      CALL GDSSUR(P(1,5),P(1,1),P(1,2),P(1,6))
      CALL GDSSUR(P(1,8),P(1,4),P(1,3),P(1,7))
      CALL GDSSUR(P(1,5),P(1,1),P(1,4),P(1,8))
      CALL GDSSUR(P(1,6),P(1,2),P(1,3),P(1,7))
C
C             Intersect cut-plane with all the surfaces of the shape
C             and draw the resulting lines
C
      CALL GDCUT
C
      GO TO 9999
C
 1500 CONTINUE
C
C             TRAP or PARA
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
C             Store all the surfaces (back,front,top,bottom,right,left)
C
      CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
      CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
      CALL GDSSUR(P(1,5),P(1,1),P(1,2),P(1,6))
      CALL GDSSUR(P(1,8),P(1,4),P(1,3),P(1,7))
      CALL GDSSUR(P(1,5),P(1,1),P(1,4),P(1,8))
      CALL GDSSUR(P(1,6),P(1,2),P(1,3),P(1,7))
C
C             Intersect cut-plane with all the surfaces of the shape
C             and draw the resulting lines
C
      CALL GDCUT
C
      GO TO 9999
C
 2000 CONTINUE
C
C             Cylindric shapes: TUBE,CONE
C
      DO 2100 I=1,40
C
        P(1,1)=RMAX1*GCOS(I)
        P(2,1)=RMAX1*GSIN(I)
        P(3,1)=Z1
        P(1,2)=RMAX1*GCOS(I+1)
        P(2,2)=RMAX1*GSIN(I+1)
        P(3,2)=Z1
        P(1,3)=RMIN1*GCOS(I+1)
        P(2,3)=RMIN1*GSIN(I+1)
        P(3,3)=Z1
        P(1,4)=RMIN1*GCOS(I)
        P(2,4)=RMIN1*GSIN(I)
        P(3,4)=Z1
        P(1,5)=RMAX2*GCOS(I)
        P(2,5)=RMAX2*GSIN(I)
        P(3,5)=Z2
        P(1,6)=RMAX2*GCOS(I+1)
        P(2,6)=RMAX2*GSIN(I+1)
        P(3,6)=Z2
        P(1,7)=RMIN2*GCOS(I+1)
        P(2,7)=RMIN2*GSIN(I+1)
        P(3,7)=Z2
        P(1,8)=RMIN2*GCOS(I)
        P(2,8)=RMIN2*GSIN(I)
        P(3,8)=Z2
C
C             Store top,bottom,back,front surfaces
C
        CALL GDSSUR(P(1,1),P(1,2),P(1,6),P(1,5))
        CALL GDSSUR(P(1,4),P(1,3),P(1,7),P(1,8))
        CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
        CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
C
 2100 CONTINUE
C
C             Intersect cut-plane with all the surfaces of the shape
C             and draw the resulting lines
C
      CALL GDCUT
C
      GO TO 9999
C
 2500 CONTINUE
C
C             Segmented cylindric shapes: TUBS,CONS
C
      CALL GDSARC(RMAX1,PHIMIN,PHIMAX,XMAX1,NP)
      CALL GDSARC(RMIN1,PHIMIN,PHIMAX,XMIN1,NP)
      CALL GDSARC(RMAX2,PHIMIN,PHIMAX,XMAX2,NP)
      CALL GDSARC(RMIN2,PHIMIN,PHIMAX,XMIN2,NP)
C
      NP1=NP-1
C
      DO 2510 I=1,NP1
C
        P(1,1)=XMAX1(1,I)
        P(2,1)=XMAX1(2,I)
        P(3,1)=Z1
        P(1,2)=XMAX1(1,I+1)
        P(2,2)=XMAX1(2,I+1)
        P(3,2)=Z1
        P(1,3)=XMIN1(1,I+1)
        P(2,3)=XMIN1(2,I+1)
        P(3,3)=Z1
        P(1,4)=XMIN1(1,I)
        P(2,4)=XMIN1(2,I)
        P(3,4)=Z1
        P(1,5)=XMAX2(1,I)
        P(2,5)=XMAX2(2,I)
        P(3,5)=Z2
        P(1,6)=XMAX2(1,I+1)
        P(2,6)=XMAX2(2,I+1)
        P(3,6)=Z2
        P(1,7)=XMIN2(1,I+1)
        P(2,7)=XMIN2(2,I+1)
        P(3,7)=Z2
        P(1,8)=XMIN2(1,I)
        P(2,8)=XMIN2(2,I)
        P(3,8)=Z2
C
C             Store top,bottom,back,front surfaces
C
        CALL GDSSUR(P(1,1),P(1,2),P(1,6),P(1,5))
        CALL GDSSUR(P(1,4),P(1,3),P(1,7),P(1,8))
        CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
        CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
C
C             Store right surface
C
        IF (I.NE.1) GO TO 2505
        CALL GDSSUR(P(1,5),P(1,1),P(1,4),P(1,8))
C
C             Store left surface
C
 2505   CONTINUE
        IF (I.NE.NP1) GO TO 2510
        CALL GDSSUR(P(1,6),P(1,2),P(1,3),P(1,7))
C
 2510 CONTINUE
C
C             Intersect cut-plane with all the surfaces of the shape
C             and draw the resulting lines
C
      CALL GDCUT
C
      GO TO 9999
C
 3000 CONTINUE
C
C             Normal view or X-view or Y-view or Z-view for SPHE
C
C             Cut not implemented
C
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
      NZ1=NZ-1
C
      DO 4050 IZ=1,NZ1
C
        ZI=PAR(5+(IZ-1)*3)
        R0=PAR(6+(IZ-1)*3)
        R1=PAR(7+(IZ-1)*3)
        ZZI=PAR(5+IZ*3)
        RR0=PAR(6+IZ*3)
        RR1=PAR(7+IZ*3)
C
        DO 4045 IDIV=1,NDIV
C
          PHI0=PHIMIN+(IDIV-1)*DPHI
          PHI1=PHI0+DPHI
          PHI0=PHI0*DEGRAD
          PHI1=PHI1*DEGRAD
C
          P(1,1)=R1*COS(PHI0)
          P(2,1)=R1*SIN(PHI0)
          P(3,1)=ZI
          P(1,2)=R1*COS(PHI1)
          P(2,2)=R1*SIN(PHI1)
          P(3,2)=ZI
          P(1,3)=R0*COS(PHI1)
          P(2,3)=R0*SIN(PHI1)
          P(3,3)=ZI
          P(1,4)=R0*COS(PHI0)
          P(2,4)=R0*SIN(PHI0)
          P(3,4)=ZI
          P(1,5)=RR1*COS(PHI0)
          P(2,5)=RR1*SIN(PHI0)
          P(3,5)=ZZI
          P(1,6)=RR1*COS(PHI1)
          P(2,6)=RR1*SIN(PHI1)
          P(3,6)=ZZI
          P(1,7)=RR0*COS(PHI1)
          P(2,7)=RR0*SIN(PHI1)
          P(3,7)=ZZI
          P(1,8)=RR0*COS(PHI0)
          P(2,8)=RR0*SIN(PHI0)
          P(3,8)=ZZI
C
C             Store top and bottom surfaces
C
          CALL GDSSUR(P(1,1),P(1,2),P(1,6),P(1,5))
          CALL GDSSUR(P(1,4),P(1,3),P(1,7),P(1,8))
C
C             Store back surface
C
          IF (IZ.NE.1) GO TO 4010
          CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
C
C             Store front surface
C
 4010     CONTINUE
          IF (IZ.NE.NZ1) GO TO 4020
          CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
C
C             Store right surface
C
 4020     CONTINUE
          IF (ABS(PHIMAX-PHIMIN).EQ.360.) GO TO 4045
          IF (IDIV.NE.1) GO TO 4030
          CALL GDSSUR(P(1,5),P(1,1),P(1,4),P(1,8))
C
C             Store left surface
C
 4030     CONTINUE
          IF (IDIV.NE.NDIV) GO TO 4045
          CALL GDSSUR(P(1,6),P(1,2),P(1,3),P(1,7))
C
 4045   CONTINUE
C
C             Intersect cut-plane with the surfaces of one Z section
C             and draw the resulting lines
C
        CALL GDCUT
C
 4050 CONTINUE
C
      GO TO 9999
C
 5000 CONTINUE
C
C             PCON
C
      NZ1=NZ-1
C
      DO 5050 IZ=1,NZ1
C
        ZI=PAR(4+(IZ-1)*3)
        R0=PAR(5+(IZ-1)*3)
        R1=PAR(6+(IZ-1)*3)
        ZZI=PAR(4+IZ*3)
        RR0=PAR(5+IZ*3)
        RR1=PAR(6+IZ*3)
C
        CALL GDSARC(R1,PHIMIN,PHIMAX,XMAX1,NP)
        CALL GDSARC(R0,PHIMIN,PHIMAX,XMIN1,NP)
        CALL GDSARC(RR1,PHIMIN,PHIMAX,XMAX2,NP)
        CALL GDSARC(RR0,PHIMIN,PHIMAX,XMIN2,NP)
        NDIV=NP-1
C
        DO 5005 IDIV=1,NDIV
C
          P(1,1)=XMAX1(1,IDIV)
          P(2,1)=XMAX1(2,IDIV)
          P(3,1)=ZI
          P(1,2)=XMAX1(1,IDIV+1)
          P(2,2)=XMAX1(2,IDIV+1)
          P(3,2)=ZI
          P(1,3)=XMIN1(1,IDIV+1)
          P(2,3)=XMIN1(2,IDIV+1)
          P(3,3)=ZI
          P(1,4)=XMIN1(1,IDIV)
          P(2,4)=XMIN1(2,IDIV)
          P(3,4)=ZI
          P(1,5)=XMAX2(1,IDIV)
          P(2,5)=XMAX2(2,IDIV)
          P(3,5)=ZZI
          P(1,6)=XMAX2(1,IDIV+1)
          P(2,6)=XMAX2(2,IDIV+1)
          P(3,6)=ZZI
          P(1,7)=XMIN2(1,IDIV+1)
          P(2,7)=XMIN2(2,IDIV+1)
          P(3,7)=ZZI
          P(1,8)=XMIN2(1,IDIV)
          P(2,8)=XMIN2(2,IDIV)
          P(3,8)=ZZI
C
C             Store top and bottom surfaces
C
          CALL GDSSUR(P(1,1),P(1,2),P(1,6),P(1,5))
          CALL GDSSUR(P(1,4),P(1,3),P(1,7),P(1,8))
C
C             Store back surface
C
          IF (IZ.NE.1) GO TO 5001
          CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
C
C             Store front surface
C
 5001     CONTINUE
          IF (IZ.NE.NZ1) GO TO 5002
          CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
C
C             Store right surface
C
 5002     CONTINUE
          IF (ABS(PHIMAX-PHIMIN).EQ.360.) GO TO 5005
          IF (IDIV.NE.1) GO TO 5003
          CALL GDSSUR(P(1,5),P(1,1),P(1,4),P(1,8))
C
C             Store left surface
C
 5003     CONTINUE
          IF (IDIV.NE.NDIV) GO TO 5005
          CALL GDSSUR(P(1,6),P(1,2),P(1,3),P(1,7))
C
 5005   CONTINUE
C
C             Intersect cut-plane with the surfaces of one Z section
C             and draw the resulting lines
C
        CALL GDCUT
C
 5050 CONTINUE
C
      GO TO 9999
*
 6000 CONTINUE
*
*          Cut tube shape.
*
      CPHIS = COS( PHIS*DEGRAD )
      SPHIS = SIN( PHIS*DEGRAD )
      P( 1, 1) = PAR(2)*CPHIS
      P( 2, 1) = PAR(2)*SPHIS
      P( 3, 1) = -( PAR( 6)*P( 1, 1) +PAR( 7)*P( 2, 1) )/PAR( 8)
     + -PAR( 3)
      P( 1, 4) = PAR(1)*CPHIS
      P( 2, 4) = PAR(1)*SPHIS
      P( 3, 4) = -( PAR( 6)*P( 1, 4) +PAR( 7)*P( 2, 4) )/PAR( 8)
     + -PAR( 3)
      P( 1, 5) = PAR(2)*CPHIS
      P( 2, 5) = PAR(2)*SPHIS
      P( 3, 5) = -( PAR( 9)*P( 1, 5) +PAR( 10)*P( 2, 5) )/PAR( 11)
     + +PAR( 3)
      P( 1, 8) = PAR(1)*CPHIS
      P( 2, 8) = PAR(1)*SPHIS
      P( 3, 8) = -( PAR( 9)*P( 1, 8) +PAR( 10)*P( 2, 8) )/PAR( 11)
     + +PAR( 3)
*
      DO 6010 I = 1, NL
        PHIS = PHIS+DPHI
        CPHIS = COS( PHIS*DEGRAD )
        SPHIS = SIN( PHIS*DEGRAD )
        P( 1, 2) = PAR(2)*CPHIS
        P( 2, 2) = PAR(2)*SPHIS
        P( 3, 2) = -( PAR( 6)*P( 1, 1) +PAR( 7)*P( 2, 1) )/PAR( 8)
     +   -PAR( 3)
        P( 1, 3) = PAR(1)*CPHIS
        P( 2, 3) = PAR(1)*SPHIS
        P( 3, 3) = -( PAR( 6)*P( 1, 4) +PAR( 7)*P( 2, 4) )/PAR( 8)
     +   -PAR( 3)
        P( 1, 6) = PAR(2)*CPHIS
        P( 2, 6) = PAR(2)*SPHIS
        P( 3, 6) = -( PAR( 9)*P( 1, 5) +PAR( 10)*P( 2, 5) )/PAR( 11)
     +   +PAR( 3)
        P( 1, 7) = PAR(1)*CPHIS
        P( 2, 7) = PAR(1)*SPHIS
        P( 3, 7) = -( PAR( 9)*P( 1, 8) +PAR( 10)*P( 2, 8) )/PAR( 11)
     +   +PAR( 3)
*
        CALL GDSSUR( P(1,1), P(1,2), P(1,6), P(1,5) )
        CALL GDSSUR( P(1,4), P(1,3), P(1,7), P(1,8) )
        CALL GDSSUR( P(1,1), P(1,2), P(1,3), P(1,4) )
        CALL GDSSUR( P(1,5), P(1,6), P(1,7), P(1,8) )
*
        IF( ISEG .EQ. 1 ) THEN
          IF( I .EQ. 1 ) THEN
            CALL GDSSUR( P(1,5), P(1,1), P(1,4), P(1,8) )
          ELSEIF ( I .EQ. NL ) THEN
            CALL GDSSUR( P(1,6), P(1,2), P(1,3), P(1,7) )
          ENDIF
        ENDIF
*
        P( 1, 1) = P( 1, 2)
        P( 2, 1) = P( 2, 2)
        P( 3, 1) = P( 3, 2)
        P( 1, 4) = P( 1, 3)
        P( 2, 4) = P( 2, 3)
        P( 3, 4) = P( 3, 3)
        P( 1, 5) = P( 1, 6)
        P( 2, 5) = P( 2, 6)
        P( 3, 5) = P( 3, 6)
        P( 1, 8) = P( 1, 7)
        P( 2, 8) = P( 2, 7)
        P( 3, 8) = P( 3, 7)
*
 6010 CONTINUE
*
      CALL GDCUT
      GO TO 9999
*
 7000 CONTINUE
C
C             ELTU
C
      DO 7010  I=1,40
C
         P(1,1)=A*GCOS(I)
         P(2,1)=B*GSIN(I)
         P(3,1)=Z1
         P(1,2)=A*GCOS(I+1)
         P(2,2)=B*GSIN(I+1)
         P(3,2)=Z1
         P(1,3)=0.
         P(2,3)=0.
         P(3,3)=Z1
         P(1,4)=A*GCOS(I)
         P(2,4)=B*GSIN(I)
         P(3,4)=Z2
         P(1,5)=A*GCOS(I+1)
         P(2,5)=B*GSIN(I+1)
         P(3,5)=Z2
         P(1,6)=0.
         P(2,6)=0.
         P(3,6)=Z2
C
C             Store top,back,front surfaces
C
         CALL GDSSUR(P(1,1),P(1,2),P(1,5),P(1,4))
         CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,3))
         CALL GDSSUR(P(1,4),P(1,5),P(1,6),P(1,6))
C
 7010 CONTINUE
C
C             Intersect cut-plane with all the surfaces of the shape
C             and draw the resulting lines
C
      CALL GDCUT
C
      GO TO 9999
C
*
 7350 CONTINUE
C
C             Another cylindrical shape: HYPE
C
      NZSTEP = 40
      DZ = Z2 / NZSTEP * 2.
      RMA1 = SQRT(RMAX12 + Z2*Z2*TANTHS)
      RMI1 = SQRT(RMIN12 + Z2*Z2*TANTHS)
      ZZ1 = Z2
      DO 7400 J = 1, NZSTEP
         ZZ2 = Z2 - J*DZ
         Z2T = (ZZ2**2) * TANTHS
         RMA2 = SQRT(RMAX12 + Z2T)
         RMI2 = SQRT(RMIN12 + Z2T)
         DO 7390  I=1,40
C
            P(1,1)=RMA1*GCOS(I)
            P(2,1)=RMA1*GSIN(I)
            P(3,1)=ZZ1
            P(1,2)=RMA1*GCOS(I+1)
            P(2,2)=RMA1*GSIN(I+1)
            P(3,2)=ZZ1
            P(1,3)=RMI1*GCOS(I+1)
            P(2,3)=RMI1*GSIN(I+1)
            P(3,3)=ZZ1
            P(1,4)=RMI1*GCOS(I)
            P(2,4)=RMI1*GSIN(I)
            P(3,4)=ZZ1
            P(1,5)=RMA2*GCOS(I)
            P(2,5)=RMA2*GSIN(I)
            P(3,5)=ZZ2
            P(1,6)=RMA2*GCOS(I+1)
            P(2,6)=RMA2*GSIN(I+1)
            P(3,6)=ZZ2
            P(1,7)=RMI2*GCOS(I+1)
            P(2,7)=RMI2*GSIN(I+1)
            P(3,7)=ZZ2
            P(1,8)=RMI2*GCOS(I)
            P(2,8)=RMI2*GSIN(I)
            P(3,8)=ZZ2
C
C             Store top and bottom surfaces only
C
            CALL GDSSUR(P(1,1),P(1,2),P(1,6),P(1,5))
            CALL GDSSUR(P(1,4),P(1,3),P(1,7),P(1,8))
C            CALL GDSSUR(P(1,1),P(1,2),P(1,3),P(1,4))
C            CALL GDSSUR(P(1,5),P(1,6),P(1,7),P(1,8))
C
 7390    CONTINUE
C
         RMA1 = RMA2
         RMI1 = RMI2
         ZZ1 = ZZ2
C
 7400  CONTINUE
C
C             Intersect cut-plane with all the surfaces of the shape
C             and draw the resulting lines
C
      CALL GDCUT
C
*
 9999 END
