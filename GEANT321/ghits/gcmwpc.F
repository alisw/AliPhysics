*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:08  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GCMWPC (DETREP,HITREP,IOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Routine to compute one or  two digitisations produced by a   *
C.    * hit on a cylindrical MWPC.                                     *
C.    * DETREP(1) number of wires                                      *
C.    * DETREP(2) wire spacing (radians)                               *
C.    * DETREP(3) dtheta/dz along the wires                            *
C.    * DETREP(4) theta of a point on wire 1                           *
C.    * DETREP(5) z of a point on wire 1                               *
C.    * DETREP(6) gap width                                            *
C.    * HITREP(1) theta coordinate of intersection                     *
C.    * HITREP(2) z coordinate                                         *
C.    * HITREP(3) dtheta/dr                                            *
C.    * HITREP(4) dz/dr                                                *
C.    * IOUT(1)   wire number (-1-missing)                             *
C.    * IOUT(2)   cluster size                                         *
C.    * IOUT(3)   wire number of second cluster if any                 *
C.    * IOUT(4)   cluster size                                         *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *       Author    M.Hansroul  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION HITREP(4), DETREP(6), IOUT(4)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION ONE
#endif
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      IOUT(1) = - 1
      IOUT(3) = - 1
      ZS     = DETREP(2)
      IF (ZS.EQ.0.)                                   GO TO 99
C.
C.            COMPUTE WIRE NUMBER
C.
      HG     = .5 * DETREP(6)
      DT     = HG * HITREP(3)
      DZ     = HG * HITREP(4)
      T1     = HITREP(1) - DT
      T2     = HITREP(1) + DT
      Z1     = HITREP(2) - DZ
      Z2     = HITREP(2) + DZ
C
      Z0     = DETREP(5)
      DT     = DETREP(3)
      DZ     = Z1 - Z0
      T1     = T1 - DZ * DT
      DZ     = Z2 - Z0
      T2     = T2 - DZ * DT
      T0     = DETREP(4)
C
      T1     = MOD (ONE*T1,TWOPI)
      T2     = MOD (ONE*T2,TWOPI)
      T0     = MOD (ONE*T0,TWOPI)
      IF (T1.LT.0.) T1 = T1 + TWOPI
      IF (T2.LT.0.) T2 = T2 + TWOPI
      IF (T0.LT.0.) T0 = T0 + TWOPI
C
      NN     = 0
      W1     = (T1 - T0) / ZS
      IF (W1.LE.0.)THEN
         NN     = NN + 1
         W1     = W1 + TWOPI / ZS
      ENDIF
C
      W2     = (T2 - T0) / ZS
      IF (W2.LE.0.)THEN
         NN     = NN + 1
         W2     = W2 + TWOPI / ZS
      ENDIF
C
      W1     = W1 + 1.5
      W2     = W2 + 1.5
      IW1    = W1
      IW2    = W2
      NWIR   = DETREP(1)
      IF (NN.NE.1)THEN
C.
C.            1 CLUSTER
C.
         IWMIN  = MIN (IW1,IW2)
         IWMAX  = MAX (IW1,IW2)
         IF (IWMIN.GT.NWIR)                              GO TO 99
C
         IWMAX  = MIN (NWIR,IWMAX)
C
         IOUT(1) = IWMIN
         IOUT(2) = IWMAX - IWMIN + 1
      ELSE
C.
C.            SPECIAL CASE: SIGNAL ON WIRE 1
C.            AND ON WIRE 'NWIRES' --> 2 CLUSTERS
C.
         IW     = MIN (IW1,IW2)
         IW     = MIN (IW,NWIR)
         IOUT(1) = 1
         IOUT(2) = IW
         IW     = MAX (IW1,IW2)
         IF (IW.GT.NWIR)                                 GO TO 99
         IOUT(3) = IW
         IOUT(4) = NWIR - IW + 1
      ENDIF
C
  99  RETURN
      END
