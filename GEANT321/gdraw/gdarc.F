*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GDARC(R,Z0,PHI1,PHI2)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw an arc, from PHI1 to PHI2,                          *
C.    *       of the circle of radius R centered on Z-axis,            *
C.    *       at a distance of Z0 from the origin,                     *
C.    *       on a plane parallel to X-Y plane                         *
C.    *                                                                *
C.    *    ==>Called by : GDRAWS                                       *
C.    *       Authors : A.McPherson, P.Zanarini   *********            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcdraw.inc"
      DIMENSION X(3,45),U(45),V(45)
C.
C.    ------------------------------------------------------------------
C.
      PHIMIN=PHI1
      PHIMAX=PHI2
C
      NP=1
      PH=PHIMIN*DEGRAD
      X(1,1)=R*COS(PH)
      X(2,1)=R*SIN(PH)
      X(3,1)=Z0
C
      IF(ABS(PHIMAX-PHIMIN).LE.9.) GO TO 20
C
      I1=PHIMIN/9.0+2
      I2=PHIMAX/9.0+1
      NAS=I2-I1+1
      IF(NAS.LE.0) NAS=NAS+40
      IC=I1
C
      DO 10 IAS=1,NAS
         IF(IC.GT.40) IC=IC-40
         NP=NP+1
         ICC=IC
         IF (ICC.LE.0) ICC=40+ICC
         X(1,NP)=R*GCOS(ICC)
         X(2,NP)=R*GSIN(ICC)
         X(3,NP)=Z0
         IC=IC+1
   10 CONTINUE
C
   20 CONTINUE
C
      PH=PHIMAX*DEGRAD
      NP=NP+1
      X(1,NP)=R*COS(PH)
      X(2,NP)=R*SIN(PH)
      X(3,NP)=Z0
C
      CALL GDFR3D(X,NP,U,V)
C
CM             Short circuit in case of GMR
C
      IF ( IGMR .EQ. 0 ) CALL GDRAWV(U,V,NP)
CM
      END
