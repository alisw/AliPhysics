*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GPMWPC (DETREP,HITREP,IOUT)
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
C     *           (6)=GAP WIDTH                                        *
C     *                                                                *
C     *     HITREP(1)=X COORDINATE OF INTERSECTION                     *
C     *           (2)=Y COORDINATE                                     *
C     *           (3)=DX/DZ                                            *
C     *           (4)=DY/DZ                                            *
C     *                                                                *
C     *     IOUT(1)=WIRE NUMBER                                        *
C     *         (2)=CLUSTER SIZE                                       *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *       Authors    F.Carena, M.Hansroul  *********               *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION HITREP(4), DETREP(6), IOUT(2)
C.
C.    ------------------------------------------------------------------
C.
      IOUT(1) = - 1
      SP     = DETREP(2)
      IF (SP.EQ.0.)                                   GO TO 99
C
      HG     = .5 * DETREP(6)
      DX     = HG * HITREP(3)
      DY     = HG * HITREP(4)
      X1     = HITREP(1) - DX
      X2     = HITREP(1) + DX
      Y1     = HITREP(2) - DY
      Y2     = HITREP(2) + DY
      SA     = DETREP(3)
      CA     = DETREP(4)
      U1     = X1 * CA + Y1 * SA
      U2     = X2 * CA + Y2 * SA
      OR     = DETREP(5)
      W1     = (U1 - OR) / SP + 1.5
      W2     = (U2 - OR) / SP + 1.5
C
      IW1    = W1
      IW2    = W2
      IWMIN  = MIN (IW1,IW2)
      IWMAX  = MAX (IW1,IW2)
      NWIR   = DETREP(1)
      IF (IWMIN.GT.NWIR)                              GO TO 99
      IF (IWMAX.LE.0)                                 GO TO 99
C
      IWMIN  = MAX (1,IWMIN)
      IWMAX  = MIN (NWIR,IWMAX)
C
      IOUT(1) = IWMIN
      IOUT(2) = IWMAX - IWMIN + 1
C
  99  RETURN
      END
