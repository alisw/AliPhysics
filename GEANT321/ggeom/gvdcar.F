*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 10/10/94  20.01.58  by  S.Giani
*-- Author :
      SUBROUTINE GVDCAR(IAXIS,ISH,IROT,PARS,CL,CH,IERR)
C.
C.    *****************************************************************
C.    *                                                               *
C.    *    ROUTINE TO FIND THE LIMITS ALONG AXIS IAXIS IN CARTESIAN   *
C.    *    COORDINATES FOR VOLUME OF SHAPE ISH ROTATED BY THE         *
C.    *    ROTATION MATRIX IROT. THE SHAPE HAS NPAR PARAMETERS IN     *
C.    *    THE ARRAY PARS. THE LOWER LIMIT IS RETURNED IN CL, THE     *
C.    *    HIGHER IN CH. IF THE CALCULATION CANNOT BE MADE IERR IS    *
C.    *    SET TO 1 OTHERWISE IT IS SET TO 0.                         *
C.    *                                                               *
C.    *    ==>Called by : GVDLIM                                      *
C.    *         Author  S.Giani  ********                             *
C.    *                                                               *
C.    *****************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
      DIMENSION PARS(50),X(3),XT(3)
C.
C.          ---------------------------------------------------
C.
      IERR=1
      IF (ISH.GT.4.AND.ISH.NE.10.AND.ISH.NE.28) GO TO 40
C
C           CUBOIDS, TRAPEZOIDS, PARALLELEPIPEDS.
C
C
      IERR=0
      CL=0
      CH=0
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
C          UPDATE LIMITS IF NECESSARY.
C
      IF(XT(IAXIS).LT.CL) CL=XT(IAXIS)
      IF(XT(IAXIS).GT.CH) CH=XT(IAXIS)
C
   30 CONTINUE
C
      GO TO 999
C
   40 CONTINUE
      IF(ISH.EQ.9) GO TO 90
C
C              TUBES , CONES, POLYGONS, POLYCONES.
C              AND CUT TUBES.
C
      MYFLAG=0
      IF((ISH.EQ.11.OR.ISH.EQ.12).AND.(IAXIS.LT.3))THEN
        MYFLAG=1
      ENDIF
      X(1)=0.0
      X(2)=0.0
      X(3)=1.0
      JROT=LQ(JROTM-IROT)
      XT(1)=X(1)
      XT(2)=X(2)
      XT(3)=X(3)
      IF(IROT.NE.0) CALL GINROT(X,Q(JROT+1),XT)
C
C          XT IS Z AXIS ROTATED.
C
      IF(MYFLAG.EQ.0)THEN
       IF(ABS(XT(IAXIS)).LT.0.99) GO TO 50
      ELSE
       IF(ABS(XT(3)).LT.0.99) GO TO 50
      ENDIF
      IF(ISH.EQ.11)GO TO 45
      IF(ISH.EQ.12)GO TO 46
C
C           PARALLEL.
C
      IP=3
      IF(ISH.GT.6.AND.ISH.NE.NSCTUB.AND.ISH.NE.13.AND.ISH.NE.14) IP=1
      CL=-PARS(IP)
      CH=PARS(IP)
      IERR=0
C
      GO TO 999
  45  IF(MYFLAG.EQ.0)THEN
       NZLAST=PARS(4)
       IZLAST=2+3*NZLAST
       CL=PARS(5)
       GO TO 49
      ELSEIF(MYFLAG.EQ.1)THEN
       NZLAST=PARS(4)
       IZLAST=2+3*NZLAST
       TMPRAD=0.
       DO 145 I=7,IZLAST+2,3
         IF(PARS(I).GT.TMPRAD)TMPRAD=PARS(I)
 145   CONTINUE
       PHIMIN=PARS(1)
       PHIMAX=PHIMIN+PARS(2)
       AANG=ABS(PHIMAX-PHIMIN)
       NANG=PARS(3)
       AATMAX=NANG*360./AANG
       LATMAX=AATMAX
       ALA=AATMAX-LATMAX
       IF(ALA.GT..5)LATMAX=LATMAX+1
       AFINV=1./COS(PI/LATMAX)
       FINV=ABS(AFINV)
       R=TMPRAD*FINV
       CL=-R
       CH= R
       IERR=0
       GOTO 999
      ENDIF
C
  46  IF(MYFLAG.EQ.0)THEN
       NZLAST=PARS(3)
       IZLAST=1+3*NZLAST
       CL=PARS(4)
      ELSEIF(MYFLAG.EQ.1)THEN
       NZLAST=PARS(3)
       IZLAST=1+3*NZLAST
       TMPRAD=0.
       DO 146 I=6,IZLAST+2,3
         IF(PARS(I).GT.TMPRAD)TMPRAD=PARS(I)
 146   CONTINUE
       CL=-TMPRAD
       CH= TMPRAD
       IERR=0
       GOTO 999
      ENDIF
C
  49  CH=PARS(IZLAST)
      IF ( ABS(XT(IAXIS)-X(IAXIS)) .GT.1.) THEN
         TEMP = CL
         CL = -CH
         CH = -TEMP
      ENDIF
      IERR=0
      GO TO 999
C
   50 CONTINUE
**
      IF(ISH.EQ.13) THEN
         CL=-PARS(IAXIS)
         CH=PARS(IAXIS)
         IERR=0
         GOTO 999
      ENDIF
**
      IF(ISH.EQ.14) THEN
C     for hyperboloid, use escribed cylinder
         CH = SQRT(PARS(2)**2+(PARS(3)*TAN(PARS(4)*DEGRAD))**2)
         CL = -CH
         IERR=0
         GOTO 999
      ENDIF
**
      IF(ISH.GT.10.AND.ISH.NE.NSCTUB)GO TO 999
      IF(ABS(XT(IAXIS)).GT.0.01) GO TO 70
C
C         Z AXIS PERPENDICULAR TO IAXIS. ASSUME COMPLETE TUBE OR
C         CONE (I.E. IGNORE PHI SEGMENTATION).
C
      IF(ISH.GT.6.AND.ISH.NE.NSCTUB) GO TO 60
C
      CL=-PARS(2)
      CH=PARS(2)
      IERR=0
      IF(ISH.EQ.6)THEN
        RMIN=PARS(1)
        RMAX=PARS(2)
        IF(IROT.NE.0)THEN
         IF(Q(JROT+15).EQ.0.)THEN
           PHI1=(PARS(4)+Q(JROT+12))*DEGRAD
           PHI2=(PARS(5)+Q(JROT+12))*DEGRAD
         ELSEIF(Q(JROT+15).EQ.180.)THEN
           PHI1=(PARS(4)+Q(JROT+12)-(PARS(5)-PARS(4)))*DEGRAD
           PHI2=(PARS(5)+Q(JROT+12)-(PARS(5)-PARS(4)))*DEGRAD
         ELSE
           GOTO 999
         ENDIF
        ELSE
         PHI1=PARS(4)*DEGRAD
         PHI2=PARS(5)*DEGRAD
        ENDIF
        IF(IAXIS.EQ.1)THEN
         IF(PHI1.GE.0..AND.PHI2.LE.PI)THEN
          XMIN1=RMIN*COS(PHI2)
          XMIN2=RMAX*COS(PHI2)
          CL=MIN(XMIN1,XMIN2)
          XMAX1=RMIN*COS(PHI1)
          XMAX2=RMAX*COS(PHI1)
          CH=MAX(XMAX1,XMAX2)
         ELSEIF(PHI1.GE.PI.AND.PHI2.LE.TWOPI.OR.
     +          PHI1.GE.-PI.AND.PHI2.LE.0.)THEN
          XMIN1=RMIN*COS(PHI1)
          XMIN2=RMAX*COS(PHI1)
          CL=MIN(XMIN1,XMIN2)
          XMAX1=RMIN*COS(PHI2)
          XMAX2=RMAX*COS(PHI2)
          CH=MAX(XMAX1,XMAX2)
         ELSEIF(PHI1.LT.0..AND.PHI2.GT.0..AND.
     +         (PHI2-PHI1).LE.PI)THEN
          XMIN1=RMIN*COS(PHI2)
          XMIN2=RMIN*COS(PHI1)
          CL1=MIN(XMIN1,XMIN2)
          XMIN3=RMAX*COS(PHI2)
          XMIN4=RMAX*COS(PHI1)
          CL2=MIN(XMIN3,XMIN4)
          CL=MIN(CL1,CL2)
          CH=RMAX
         ELSEIF(PHI1.LT.PI.AND.PHI2.GT.PI.AND.
     +         (PHI2-PHI1).LE.PI)THEN
          CL=-RMAX
          XMAX1=RMIN*COS(PHI2)
          XMAX2=RMIN*COS(PHI1)
          CH1=MAX(XMAX1,XMAX2)
          XMAX3=RMAX*COS(PHI2)
          XMAX4=RMAX*COS(PHI1)
          CH2=MAX(XMAX3,XMAX4)
          CH=MAX(CH1,CH2)
         ENDIF
        ELSEIF(IAXIS.EQ.2)THEN
         IF(PHI1.GE.(-PI*.5).AND.PHI2.LE.(PI*.5))THEN
          YMIN1=RMIN*SIN(PHI1)
          YMIN2=RMAX*SIN(PHI1)
          CL=MIN(YMIN1,YMIN2)
          YMAX1=RMIN*SIN(PHI2)
          YMAX2=RMAX*SIN(PHI2)
          CH=MAX(YMAX1,YMAX2)
         ELSEIF(PHI1.GE.(PI*.5).AND.PHI2.LE.(PI*3*.5))THEN
          YMIN1=RMIN*SIN(PHI2)
          YMIN2=RMAX*SIN(PHI2)
          CL=MIN(YMIN1,YMIN2)
          YMAX1=RMIN*SIN(PHI1)
          YMAX2=RMAX*SIN(PHI1)
          CH=MAX(YMAX1,YMAX2)
         ELSEIF(PHI1.LT.(PI*.5).AND.PHI2.GT.(PI*.5).AND.
     +         (PHI2-PHI1).LE.PI)THEN
          YMIN1=RMIN*SIN(PHI2)
          YMIN2=RMIN*SIN(PHI1)
          CL1=MIN(YMIN1,YMIN2)
          YMIN3=RMAX*SIN(PHI2)
          YMIN4=RMAX*SIN(PHI1)
          CL2=MIN(YMIN3,YMIN4)
          CL=MIN(CL1,CL2)
          CH=RMAX
         ELSEIF(((PHI1.LT.(PI*3*.5).AND.PHI2.GT.(PI*3*.5)).OR.
     +          (PHI1.LT.-(PI*.5).AND.PHI2.GT.-(PI*.5)))
     +          .AND.(PHI2-PHI1).LE.PI)THEN
          CL=-RMAX
          YMAX1=RMIN*SIN(PHI2)
          YMAX2=RMIN*SIN(PHI1)
          CH1=MAX(YMAX1,YMAX2)
          YMAX3=RMAX*SIN(PHI2)
          YMAX4=RMAX*SIN(PHI1)
          CH2=MAX(YMAX3,YMAX4)
          CH=MAX(CH1,CH2)
         ENDIF
        ENDIF
      ENDIF
C
      GO TO 999
C
   60 CONTINUE
C
      RM=PARS(3)
      IF(PARS(5).GT.PARS(3)) RM=PARS(5)
C
      CL=-RM
      CH=RM
      IERR=0
C
      GO TO 999
C
   70 CONTINUE
C
C           ARBITRARY ROTATION.
C
      DZ=PARS(3)
      RM=PARS(2)
      IF(ISH.EQ.13) THEN
**
**       approxime to a cylinder whit radius
**       equal to the ellipse major axis
**
         IF(PARS(1).GT.RM) RM=PARS(1)
         GOTO 80
      ENDIF
**
      IF(ISH.EQ.14) THEN
        RM = SQRT(PARS(2)**2+(PARS(3)*TAN(PARS(4)*DEGRAD))**2)
        GO TO 80
      ENDIF
*
      IF(ISH.EQ.NSCTUB) THEN
        S1 = (1.0-PARS(8))*(1.0+PARS(8))
        IF( S1 .GT. 0.0) S1 = SQRT(S1)
        S2 = (1.0-PARS(11))*(1.0+PARS(11))
        IF( S2 .GT. 0.0) S2 = SQRT(S2)
        IF( S2 .GT. S1 ) S1 = S2
        DZ = DZ+RM*S1
      ENDIF
      IF(ISH.LE.6) GO TO 80
C
      DZ=PARS(1)
      RM=PARS(3)
      IF(PARS(5).GT.RM) RM=PARS(5)
C
   80 CONTINUE
C
      COST=ABS(XT(IAXIS))
      SINT=(1+COST)*(1-COST)
      IF(SINT.GT.0.0) SINT=SQRT(SINT)
C
      CH=COST*DZ+SINT*RM
      CL=-CH
      IERR=0
C
      GO TO 999
   90 CONTINUE
C
C           SPHERE - ASSUME COMPLETE SPHERE, TAKE OUTER RADIUS.
C
      IERR=0
      CL=-PARS(2)
      CH=PARS(2)
C
  999 CONTINUE
      END
