*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:35  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GVROT(DCOSIN,PART)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Rotates vector from one reference system into another         *
C.    *  defined by the direction cosines contained in the vector      *
C.    *  DCOSIN. The input vector PART is overwritten.                 *
C.    *                                                                *
C.    *    ==>Called by :                                              *
C.    *       Author    F.Carminati **********                         *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION PART(*),DCOSIN(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION DUX,DUY,DUZ,SINTH2,ONE,COSTH,SINTH,COSPH,SINPH
      DOUBLE PRECISION DNORM
#endif
      LOGICAL ROTATE
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      DUX      = DCOSIN(1)
      DUY      = DCOSIN(2)
      DUZ      = DCOSIN(3)
      ROTATE   = .TRUE.
      IF(ABS(DUZ).GE.0.85) THEN
         SINTH2=DUX**2+DUY**2
         IF(SINTH2.GT.0.) THEN
            COSTH  = SIGN(ONE,DUZ)*SQRT(ONE-SINTH2)
            SINTH  = SQRT(SINTH2)
            COSPH  = DUX/SINTH
            SINPH  = DUY/SINTH
         ELSEIF(DUZ.GT.0.)THEN
            ROTATE = .FALSE.
         ELSE
            COSTH = -1.
            SINTH =  0.
            COSPH =  1.
            SINPH =  0.
         ENDIF
      ELSE
         COSTH  = DUZ
         SINTH  = SQRT((ONE+DUZ)*(ONE-DUZ))
         DNORM  = ONE/SQRT(DUX**2+DUY**2)
         COSPH  = DUX*DNORM
         SINPH  = DUY*DNORM
      ENDIF
*
      IF(ROTATE) THEN
         P1=PART(1)
         P2=PART(2)
         P3=PART(3)
         PART(1)= P1*COSTH*COSPH - P2*SINPH + P3*SINTH*COSPH
         PART(2)= P1*COSTH*SINPH + P2*COSPH + P3*SINTH*SINPH
         PART(3)=-P1*SINTH                  + P3*COSTH
      ENDIF
C
      END
