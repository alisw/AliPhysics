*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 06/07/94  18.28.55  by  S.Giani
*-- Author :
      SUBROUTINE GFANG(P,COSTH,SINTH,COSPH,SINPH,ROTATE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Finds theta-phi angles for a particle with director cosines P *
C.    *  The logical flag ROTATE indicates whether the the director    *
C.    *  cosines are along the Z axis (.FALSE.) or not (.TRUE.).       *
C.    *                                                                *
C.    *    ==>Called by :                                              *
C.    *       Author    F.Carminati  *********                         *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION P(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION DUX, DUY, DUZ, DSITH2, ONE, DNORM
#endif
      LOGICAL ROTATE
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      DUX      = P(1)
      DUY      = P(2)
      DUZ      = P(3)
      ROTATE   = .TRUE.
      IF(ABS(DUZ).GE.0.85) THEN
         DSITH2=DUX**2+DUY**2
         IF(DSITH2.GT.0.) THEN
            COSTH  = SIGN(ONE,DUZ)*SQRT(ONE-DSITH2)
            DSITH  = SQRT(DSITH2)
            SINTH  = DSITH
            COSPH  = DUX/DSITH
            SINPH  = DUY/DSITH
         ELSEIF(DUZ.GT.0.)THEN
            ROTATE = .FALSE.
            COSTH =  1.
            SINTH =  0.
            COSPH =  1.
            SINPH =  0.
         ELSE
            COSTH = -1.
            SINTH =  0.
            COSPH =  1.
            SINPH =  0.
         ENDIF
      ELSE
         COSTH  = DUZ
         DSITH  = SQRT((ONE+DUZ)*(ONE-DUZ))
         SINTH  = DSITH
         DNORM  = ONE/SQRT(DUX**2+DUY**2)
         COSPH  = DUX*DNORM
         SINPH  = DUY*DNORM
      ENDIF
*
      END
