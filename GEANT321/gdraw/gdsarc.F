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
      SUBROUTINE GDSARC(R,PHI1,PHI2,X,NP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Computes and store in X(1..2,1..NP) arc points           *
C.    *       (only X and Y coordinates) for GDRWSC                    *
C.    *       of the circle of radius R centered on Z-axis,            *
C.    *       from PHI1 to PHI2,                                       *
C.    *       on a plane parallel to X-Y plane                         *
C.    *                                                                *
C.    *    ==>Called by : GDRWSC                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
      DIMENSION X(2,42)
      SAVE MAXNP
      DATA MAXNP/42/
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
C
      IF(ABS(PHIMAX-PHIMIN).LE.9.) GO TO 20
C
      I1=PHIMIN/9.0+2
      I2=PHIMAX/9.0+1
      NAS=I2-I1+1
      IF(NAS.LE.0) NAS=NAS+40
      IC=I1
C
      IF (NAS.GT.MAXNP-2) GO TO 30
C
      DO 10 IAS=1,NAS
         IF(IC.GT.40) IC=IC-40
         NP=NP+1
         ICC=IC
         IF (ICC.LE.0) ICC=40+ICC
         X(1,NP)=R*GCOS(ICC)
         X(2,NP)=R*GSIN(ICC)
         IC=IC+1
   10 CONTINUE
C
   20 CONTINUE
C
      PH=PHIMAX*DEGRAD
      NP=NP+1
      X(1,NP)=R*COS(PH)
      X(2,NP)=R*SIN(PH)
      GO TO 999
C
   30 CONTINUE
      NP=NAS+2
      WRITE (CHMAIL,1000) NP
      CALL GMAIL(0,0)
C
 1000 FORMAT (' GDSARC: INTERNAL ERROR !!!   -   NP =',I5)
  999 RETURN
      END
