*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GPDRIF (DETREP,HITREP,IOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C     *     DETREP(1)=NUMBER OF WIRES                                  *
C     *           (2)=WIRE SPACING                                     *
C     *           (3)=SIN (ALPHA) (ALPHA=ANGLE OF THE NORMAL           *
C     *                            TO THE WIRE WITH RESPECT TO         *
C     *                            AXIS I)                             *
C     *           (4)=COS (ALPHA)                                      *
C     *           (5)=DISTANCE OF WIRE 1 FROM THE ORIGIN               *
C     *           (6)=DRIFT VELOCITY (CM/NSEC)                         *
C     *                                                                *
C     *     HITREP(1)=X COORDINATE OF INTERSECTION                     *
C     *           (2)=Y COORDINATE                                     *
C     *                                                                *
C     *     IOUT(1)=WIRE NUMBER                                        *
C     *         (2)=DRIFT TIME (SIGNED TO AVOID RIGHT/LEFT AMBIGUITY)  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *       Authors    F.Carena, M.Hansroul  *********               *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION HITREP(2), DETREP(6), IOUT(2)
C.
C.    ------------------------------------------------------------------
C.
      IOUT(1) = - 1
      SP     = DETREP(2)
      DV     = DETREP(6)
      IF (SP.EQ.0.)                                   GO TO 99
      IF (DV.EQ.0.)                                   GO TO 99
C
      X      = HITREP(1)
      Y      = HITREP(2)
      SA     = DETREP(3)
      CA     = DETREP(4)
      U      = X * CA + Y * SA
      OR     = DETREP(5)
      W      = (U - OR) / SP + 1.5
C
      IW     = W
      NWIR   = DETREP(1)
      IF (IW.GT.NWIR)                                 GO TO 99
      IF (IW.LE.0)                                    GO TO 99
C
      DIST   = U - OR - (IW - 1) * SP
      TDRIF  = DIST / DV
C
      IOUT(1) = IW
      IOUT(2) = TDRIF
C
  99  RETURN
      END
