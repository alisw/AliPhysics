*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GKXYZ(EPSILO)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints the kinematics of the track picked (using the     *
C.    *       graphics cursor) in case that the space points of the    *
C.    *       track had been stored in the bank JXYZ.                  *
C.    *       The picking is repeated until the character typed is     *
C.    *       'Q' or 'q'.                                              *
C.    *       EPSILO is the delta angle used for pick, in the sense    *
C.    *       that tracks to be considered as possible candidates      *
C.    *       must start with a direction in the angular range         *
C.    *       between ANG0-EPSILO and ANG0+EPSILO (where ANG0 is the   *
C.    *       angle relative to the cursor position).                  *
C.    *       If EPSILO is zero there is no delta angle concept and    *
C.    *       over all the track points the one nearest to the pick    *
C.    *       point is taken.                                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
      DIMENSION U(50),V(50)
      DIMENSION UP(5),VP(5)
C
      DIMENSION XZERO(3)
C
      SAVE XZERO
      DATA XZERO/3*0./
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('OFF ',IOFF,4,4)
C
C            Set IOBJ to TRACK
C
      IOBJ=2
C
   10 CONTINUE
      CALL GDCURS(U0,V0,KCHAR)
      IF (KCHAR.EQ.0) GO TO 130
C
      IF(JXYZ.LE.0)GO TO 130
      NXYZ=IQ(JXYZ-2)
      N1  =1
      N2  =NXYZ
C
      CALL GDFR3D(XZERO,1,UOR,VOR)
      VATAN=V0-VOR
      UATAN=U0-UOR
      IF (VATAN.EQ.0..AND.UATAN.EQ.0.) THEN
         WRITE (CHMAIL,2200)
         CALL GMAIL(0,0)
         GO TO 130
      ENDIF
      ANG0=ATAN2(VATAN,UATAN)
      MAXWS=10000
      IF (EPSILO.EQ.0) MAXWS=MAXWS*2
 
C
   20 CONTINUE
      CALL GWORK(MAXWS)
      ANGEPS=ABS(EPSILO)
      EPSIL=1.
      EPS=MAX(PLTRNX,PLTRNY)*EPSIL/20.
      IBIAS=0
C
   30 CONTINUE
C
      DO 90 I=N1,N2
C
C             Track is not (yet) a canditate
C
         IFL=0
         ILOOP=0
C
         JX=LQ(JXYZ-I)
         IF(JX.LE.0)GO TO 90
         IPOINT=JX+3
   40    IF(IPOINT.GT.JX+IQ(JX-1))GO TO 90
         NPOINT=Q(IPOINT)
         IF(NPOINT.LE.0)GO TO 90
         IPF=1
   50    IPL=IPF+49
         IF(IPL.GT.NPOINT)IPL=NPOINT
         NP=IPL-IPF+1
         IF(NP.LE.1)GO TO 80
C
C            Is THRZ option set on (R-Z projection)
C            and the cut lateral (vertical or horizontal) ?
C
         IF (ITHRZ.NE.IOFF) THEN
            IF (ICUT.EQ.1.OR.ICUT.EQ.2) ITR3D=-I
         ELSE
            ITR3D=0
         ENDIF
C
         IF (EPSILO.EQ.0) GO TO 60
C
         IF (ILOOP.EQ.0) THEN
            ILOOP=1
C
            CALL GDFR3D(Q(IPOINT+3*IPF-1),2,U,V)
            CALL GDFR3D(Q(IPOINT+3*NPOINT-1),1,ULAST,VLAST)
C
            IF (U0.GE.UOR-EPS.AND.U0.LE.UOR+EPS.AND. V0.GE.VOR-EPS.AND.
     +      V0.LE.VOR+EPS) THEN
C
C             Pick point too near the origin
C
               IFL=1
               GO TO 60
            ELSE IF (U(1).LT.UOR-EPS.OR.U(1).GT.UOR+EPS.OR. V(1).LT.VOR-
     +      EPS.OR.V(1).GT.VOR+EPS) THEN
C
C             First track point does not start near the origin
C
               IFL=1
               GO TO 60
            ELSE
               VATAN=V(2)-V(1)
               UATAN=U(2)-U(1)
               IF (VATAN.EQ.0..AND.UATAN.EQ.0.) THEN
                  IFL=1
                  GO TO 60
               ENDIF
               ANG=ATAN2(VATAN,UATAN)
               VATAN=VLAST-V(1)
               UATAN=ULAST-U(1)
               IF (VATAN.EQ.0..AND.UATAN.EQ.0.) THEN
                  IFL=1
                  GO TO 60
               ENDIF
               ANGLST=ATAN2(VATAN,UATAN)
               IF (ANG.GT.(PI-ANGEPS).OR.ANG.LT.(-PI+ANGEPS)) THEN
C
C             First track point near the +PI/-PI edge
C
                  IFL=1
                  GO TO 60
               ELSE IF (ANG.GT.ANG0-ANGEPS.AND.ANG.LT.ANG0+ANGEPS) THEN
C
C             Pick point inside the pick angle
C
                  IFL=1
                  GO TO 60
               ELSE IF (ANGLST.LT.ANG-ANGEPS.OR.ANGLST.GT.ANG+ANGEPS)
     +         THEN
C
C             Track is too curved
C
                  IFL=1
                  GO TO 60
               ENDIF
            ENDIF
C
         ENDIF
C
C             Skip completely the track ?
C
         IF (IFL.NE.1) GO TO 90
C
   60    CONTINUE
C
         IF (IBIAS+3*NP.GE.MAXWS) THEN
            MAXWS=MAXWS*2
            GO TO 20
         ENDIF
C
         CALL GDFR3D(Q(IPOINT+3*IPF-1),NP,U,V)
C
         DO 70 L=1,NP
            WS(3*L-2+IBIAS)=U(L)
            WS(3*L-1+IBIAS)=V(L)
            WS(3*L+IBIAS)=I
   70    CONTINUE
C
         IBIAS=IBIAS+3*NP
C
         IPF=IPL
         IF(IPF.LT.NPOINT)GO TO 50
   80    IPOINT=IPOINT+3.*Q(IPOINT)+2.
         GO TO 40
C
   90 CONTINUE
C
C
      IF (IBIAS.EQ.0) THEN
         ANGEPS=ANGEPS*2.
         GO TO 30
      ENDIF
C
      DISMIN=1.
      DISTPK=1000000.
      IPK=0
C
      DO 120 I=1,IBIAS,3
         UPK=WS(I)
         VPK=WS(I+1)
         DIST=(UPK-U0)*(UPK-U0)+(VPK-V0)*(VPK-V0)
         IF (DIST.LT.DISTPK) THEN
            DISTPK=DIST
            IPK=I
         ENDIF
  120 CONTINUE
C
      IF (DISTPK.GT.DISMIN) THEN
         IT=0
      ELSE
         UPK=WS(IPK)
         VPK=WS(IPK+1)
         IT=WS(IPK+2)
      ENDIF
C
      IF (IT.EQ.0) THEN
C
         WRITE (CHMAIL,2000)
         CALL GMAIL(0,0)
C
      ELSE
C
C             Draw the picked track point
C
         UDELTA=0.05*PLTRNX/20.
         VDELTA=0.05*PLTRNY/20.
         UP(1)=UPK-UDELTA
         VP(1)=VPK-VDELTA
         UP(2)=UPK+UDELTA
         VP(2)=VP(1)
         UP(3)=UP(2)
         VP(3)=VPK+VDELTA
         UP(4)=UP(1)
         VP(4)=VP(3)
         UP(5)=UP(1)
         VP(5)=VP(1)
         CALL GDRAWV(UP,VP,5)
C
C             Print kinematics of the picked track
C
         CALL GPKINE(IT)
C
      ENDIF
C
C             Restart for another pick operation
C
      GO TO 10
C
  130 CONTINUE
C
C             Reset ITR3D to 'standard projection'
C
      ITR3D=0
      IOBJ=0
      RETURN
 1900 FORMAT (' GKXYZ: Work space doubled - Now is',I6)
 2000 FORMAT (' GKXYZ: No track picked')
 2100 FORMAT (' GKXYZ: Pick angle doubled - Now is',F5.2)
 2200 FORMAT (' GKXYZ: Primary vertex: no track picked')
      END
