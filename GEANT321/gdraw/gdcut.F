*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCUT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the 4 lines obtained by the intersection            *
C.    *       of the cut plane (at point DCUT along axis ICUT) and     *
C.    *       the NSURF surfaces stored in vector SURF.                *
C.    *       DCUT, ICUT, NSURF, SURF are in common block GCDRAW.      *
C.    *                                                                *
C.    *    ==>Called by : GDRWSC, GDSSUR                               *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      DIMENSION SURF(2)
      EQUIVALENCE (WS(1),SURF(1))
      DIMENSION P(3,2),U(2),V(2),XI(2),YI(2),ZI(2),PL(3),PR(3)
      EQUIVALENCE (PL(1),XL),(PL(2),YL),(PL(3),ZL)
      EQUIVALENCE (PR(1),XR),(PR(2),YR),(PR(3),ZR)
C
      SAVE EPS,ICUT3
      DATA EPS/0.00001/
C.
C.    ------------------------------------------------------------------
C.
C     Let be L(XL,YL,ZL) and R(XR,YR,ZR) the extreme points of a line,
C     I(XI,YI,ZI) the intersection point between the line and
C     the cut-plane (for example Z=C, i.e. the plane normal to Z axis).
C     Then we have:
C
C
C                       (XR-XL)               (XR-XL)
C     XI = XL + (ZI-ZL) ------- = XL + (C-ZL) -------
C                       (ZR-ZL)               (ZR-ZL)
C
C
C                       (YR-YL)               (YR-YL)
C     YI = YL + (XI-XL) ------- = YL + (C-XL) -------
C                       (XR-XL)               (XR-XL)
C
C
C     ZI = C
C
C
C     Then we have two cases:
C
C     1.   ZL <= ZI <= ZR
C          i.e. the cut-plane intersects the line at I(XI,YI,ZI) :
C          a line is drawn between that point and next one
C
C     2.   ZL = ZR
C          i.e. the line is parallel to (or into) the cut-plane :
C          any line is drawn
C
C
      IF (NSURF.EQ.0) GO TO 80
C
      IF (ICUT.EQ.0) THEN
         ICUT3=3
         ICUT=3
      ELSE
         ICUT3=0
      ENDIF
C
      IF (ICUT.LT.1.OR.ICUT.GT.3) GO TO 80
C
      DO 70  I=1,NSURF
C
         ISTART=12*(I-1)+1
         J=0
C
         DO 40 L=1,4
C
            IL=ISTART+(L-1)*3
            IR=IL+3
            IF (L.EQ.4) IR=ISTART
            XL=SURF(IL)
            YL=SURF(IL+1)
            ZL=SURF(IL+2)
            XR=SURF(IR)
            YR=SURF(IR+1)
            ZR=SURF(IR+2)
C
            IF (ICUT3.GT.0) THEN
C
C             Transform PL and PR to have Z-axis = line of sight CTHETA,
C
               SAVTHE=GTHETA
               SAVPHI=GPHI
               GTHETA=CTHETA
               GPHI=CPHI
               NGVIEW=0
               CALL GD3D3D(PL,-1,PL,0)
               CALL GD3D3D(PR,-1,PR,0)
               GTHETA=SAVTHE
               GPHI=SAVPHI
               NGVIEW=0
C
            ENDIF
C
C             Is the line parallel to the cut-plane ?
C
            IF ((PL(ICUT).GE.PR(ICUT)-EPS).AND.(PL(ICUT).LE.PR(ICUT)+
     +      EPS)) GO TO 40
C
C             No, it isn't
C
            AMX=MAX(PL(ICUT),PR(ICUT))
            AMN=MIN(PL(ICUT),PR(ICUT))
            IF ((DCUT.LT.AMN).OR.(DCUT.GT.AMX)) GO TO 40
C
            J=J+1
C
            IF (ICUT.NE.3) GO TO 10
            XI(J)=XL+(DCUT-ZL)*(XR-XL)/(ZR-ZL)
            YI(J)=YL+(DCUT-ZL)*(YR-YL)/(ZR-ZL)
            ZI(J)=DCUT
            GO TO 30
C
   10       IF (ICUT.NE.1) GO TO 20
            XI(J)=DCUT
            YI(J)=YL+(DCUT-XL)*(YR-YL)/(XR-XL)
            ZI(J)=ZL+(DCUT-XL)*(ZR-ZL)/(XR-XL)
            GO TO 30
C
   20       CONTINUE
            XI(J)=XL+(DCUT-YL)*(XR-XL)/(YR-YL)
            YI(J)=DCUT
            ZI(J)=ZL+(DCUT-YL)*(ZR-ZL)/(YR-YL)
C
   30       IF (J.EQ.2) GO TO 50
C
   40    CONTINUE
C
         IF (J.LT.2) GO TO 70
C
   50    CONTINUE
         DO 60 IND=1,2
            P(1,IND)=XI(IND)
            P(2,IND)=YI(IND)
            P(3,IND)=ZI(IND)
   60    CONTINUE
C
         IF (ICUT3.GT.0) THEN
C
C             Transform back P
C
            SAVTHE=GTHETA
            SAVPHI=GPHI
            GTHETA=CTHETA
            GPHI=CPHI
            NGVIEW=0
            CALL GD3D3D(P(1,1),-2,P(1,1),1)
            GTHETA=SAVTHE
            GPHI=SAVPHI
            NGVIEW=0
C
         ENDIF
C
C             Points are now in general reference system
C
         CALL GDFR3D(P,-2,U,V)
C
CM             Short circuit in case of GMR
C
         IF ( IGMR .EQ. 0 ) CALL GDRAWV(U,V,2)
CM
   70 CONTINUE
C
   80 CONTINUE
C
      IF (ICUT3.GT.0) THEN
         ICUT3=0
         ICUT=0
      ENDIF
C
      NSURF=0
      ISURF=1
C
      END
