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
      SUBROUTINE GFLRAD(IAXIS,ISH,IROT,DX,PARS,CL,CH,IERR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO COMPUTE THE LIMITS IN R FOR THE SHAPE ISH        *
C.    *    DISPLACED BY THE VECTOR DX AND ROTATED BY THE MATRIX IROT.  *
C.    *    IF IAXIS = 4 THE R IS THE XY PLANE R, IF IAXIS = 5 IT IS    *
C.    *    THE 3 DINEMSIONAL SPACE R. THE SHAPE HAS NPAR PARAMETERS    *
C.    *    IN THE ARRAY PARS. THE LOWER LIMIT IS RETURNED IN CL AND    *
C.    *    THE HIGHER IN CH. IF THE CALCULATION CANNOT BE PERFORMED    *
C.    *    IERR IS SET TO 1 OTHERWISE IT IS SET TO 0.                  *
C.    *                                                                *
C.    *    ==>Called by : GFCLIM                                       *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
      DIMENSION DX(3),PARS(11),X(3),XT(3)
C.
C.           --------------------------------------------------
C.
      IERR=1
C
C            FIRST CALCULATE THE LENGTH OF THE DISPLACEMENT OF THE
C            ORIGIN.
C
      DXS=DX(1)*DX(1)+DX(2)*DX(2)
      IF(IAXIS.EQ.5) DXS=DXS+DX(3)*DX(3)
      IF(DXS.GT.0.0) DXS=SQRT(DXS)
C
      IF(ISH.GT.4.AND.ISH.NE.10.AND.ISH.NE.28) GO TO 40
C
C          CUBOIDS, TRAPEZOIDS, PARALLELEPIPEDS.
C
      CH=0.0
      CL=DXS
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
C          NOW COMPUTE RMIN = PROJECTED R ON DX AND RMAX = R
C          AND UPDATE LIMITS IF NECESSARY.
C
      R2=(XT(1)+DX(1))**2+(XT(2)+DX(2))**2
      IF(IAXIS.EQ.5) R2=R2+(XT(3)+DX(3))**2
      R=SQRT(R2)
      IF(R.GT.CH) CH=R
C
      IF(CL.LE.0.0) GO TO 30
C
      XPT=DX(1)*XT(1)+DX(2)*XT(2)
      IF(IAXIS.EQ.5) XPT=XPT+DX(3)*XT(3)
      IF(DXS.LE.1.0E-05) GO TO 30
      RMN=DXS+XPT/DXS
      IF(RMN.LT.CL) CL=RMN
C
   30 CONTINUE
C
      IF(CL.LE.0.0) CL=0.0
C
      IERR=0
      GO TO 999
C
   40 CONTINUE
      IF(ISH.GT.8.AND.ISH.NE.NSCTUB.AND.ISH.NE.13.AND.ISH.NE.14)GO TO 80
C
C             TUBES AND CONES.
C
      IP3=3
      IF(ISH.GT.6.AND.ISH.NE.NSCTUB.AND.ISH.NE.13.AND.ISH.NE.14) IP3=1
      DZ=PARS(IP3)
      R=PARS(2)
      IF(ISH.EQ.NSCTUB) THEN
        S1 = (1.0-PARS(8))*(1.0+PARS(8))
        IF( S1 .GT. 0.0) S1 = SQRT(S1)
        S2 = (1.0-PARS(11))*(1.0+PARS(11))
        IF( S2 .GT. 0.0) S2 = SQRT(S2)
        IF( S2 .GT. S1 ) S1 = S2
        DZ = DZ+R*S1
      ENDIF
**
      IF(ISH.EQ.13) THEN
**
**       APPROXIME TO A CYLINDER WHIT RADIUS
**       EQUAL TO THE ELLIPSE MAJOR AXIS
**
         RMN=0.0
         IF(PARS(1).GT.R) R=PARS(1)
         GOTO 50
      ENDIF
      RMN=PARS(1)
*
      IF(ISH.EQ.14) THEN
        R = SQRT(PARS(2)**2+(PARS(3)*TAN(PARS(4)*DEGRAD))**2)
        GO TO 50
      ENDIF
C
      IF(ISH.LE.6.OR.ISH.EQ.NSCTUB) GO TO 50
C
      R=PARS(3)
      IF(PARS(5).GT.R) R=PARS(5)
      RMN=PARS(2)
      IF(PARS(4).LT.RMN) RMN=PARS(4)
C
   50 CONTINUE
C
C          ROTATE THE LOCAL Z AXIS.
C
      X(1)=0.0
      X(2)=0.0
      X(3)=1.0
      JROT=LQ(JROTM-IROT)
      XT(1)=X(1)
      XT(2)=X(2)
      XT(3)=X(3)
      IF(IROT.NE.0) CALL GINROT(X,Q(JROT+1),XT)
C
C          COMPUTE RMIN AND RMAX ASSUMING COMPLETE TUBE HALF
C          LENGTH DZ AND RADIUS R.
C
      ST2=1.0
      IF(IAXIS.EQ.4) ST2=(1+XT(3))*(1-XT(3))
      DR=SQRT(DZ*DZ*ST2+R*R)
      CL=DXS-DR
      IF(CL.LT.0.0) CL=0.0
      CH=DXS+DR
      IF(IROT.EQ.0.AND.DXS.LT.1.0E-05) CL=RMN
      IERR=0
C
      GO TO 999
C
   80 CONTINUE
      IF(ISH.GT.9) GO TO 999
C
C           SPHERE.
C
      CL=DXS-PARS(2)
      IF(CL.LT.0.0) CL=0.0
      CH=DXS+PARS(2)
      IF(IAXIS.EQ.5.AND.DXS.LT.1.0E-05) CL=PARS(1)
      IERR=0
C
  999 CONTINUE
      END
