*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GHELIX (CHARGE, STEP, VECT, VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Performs the tracking of one step in a magnetic field         *
C.    *  The trajectory is assumed to be a helix in a constant field   *
C.    *  taken at the mid point of the step.                           *
C.    *  Parameters:                                                   *
C.    *   input                                                        *
C.    *     STEP =arc length of the step asked                         *
C.    *     VECT =input vector (position,direction cos and momentum)   *
C.    *     CHARGE=  electric charge of the particle                   *
C.    *   output                                                       *
C.    *     VOUT = same as VECT after completion of the step           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Author    M.Hansroul  *********                          *
C.    *       Modified  S.Egli, S.V.Levonian                           *
C.    *       Modified  V.Perevoztchikov
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION      VECT(7),VOUT(7)
      DIMENSION      XYZ(3),H(4),HXP(3)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6,IPP=7)
      PARAMETER      (SIXTH = 1./6.)
      PARAMETER      (EC=2.9979251E-4)
C.
C.    ------------------------------------------------------------------
C.
C       units are kgauss,centimeters,gev/c
C
      VOUT(IPP) = VECT(IPP)
      IF (CHARGE.EQ.0.)         GO TO 10
      XYZ(1)    = VECT(IX) + 0.5 * STEP * VECT(IPX)
      XYZ(2)    = VECT(IY) + 0.5 * STEP * VECT(IPY)
      XYZ(3)    = VECT(IZ) + 0.5 * STEP * VECT(IPZ)
C
      CALL GUFLD (XYZ, H)
 
      H2XY = H(1)**2 + H(2)**2
      H(4) = H(3)**2 + H2XY
      IF (H(4).LE.1.E-12)       GO TO 10
      IF (H2XY.LE.1.E-12*H(4))  THEN
         CALL GHELX3 (CHARGE*H(3), STEP, VECT, VOUT)
         GO TO 999
      ENDIF
      H(4) = SQRT(H(4))
      H(1) = H(1)/H(4)
      H(2) = H(2)/H(4)
      H(3) = H(3)/H(4)
      H(4) = H(4)*EC
*
      HXP(1) = H(2)*VECT(IPZ) - H(3)*VECT(IPY)
      HXP(2) = H(3)*VECT(IPX) - H(1)*VECT(IPZ)
      HXP(3) = H(1)*VECT(IPY) - H(2)*VECT(IPX)
 
      HP = H(1)*VECT(IPX) + H(2)*VECT(IPY) + H(3)*VECT(IPZ)
*
      RHO = -CHARGE*H(4)/VECT(IPP)
      TET = RHO * STEP
      IF (ABS(TET).GT.0.15)     THEN
         SINT = SIN(TET)
         SINTT = (SINT/TET)
         TSINT = (TET-SINT)/TET
         COS1T = 2.*(SIN(0.5*TET))**2/TET
      ELSE
         TSINT = SIXTH*TET**2
         SINTT = (1. - TSINT)
         SINT = TET*SINTT
         COS1T = 0.5*TET
      ENDIF
*
      F1 = STEP * SINTT
      F2 = STEP * COS1T
      F3 = STEP * TSINT * HP
      F4 = -TET*COS1T
      F5 = SINT
      F6 = TET * COS1T * HP
 
      VOUT(IX) = VECT(IX) + (F1*VECT(IPX) + F2*HXP(1) + F3*H(1))
      VOUT(IY) = VECT(IY) + (F1*VECT(IPY) + F2*HXP(2) + F3*H(2))
      VOUT(IZ) = VECT(IZ) + (F1*VECT(IPZ) + F2*HXP(3) + F3*H(3))
 
      VOUT(IPX) = VECT(IPX) + (F4*VECT(IPX) + F5*HXP(1) + F6*H(1))
      VOUT(IPY) = VECT(IPY) + (F4*VECT(IPY) + F5*HXP(2) + F6*H(2))
      VOUT(IPZ) = VECT(IPZ) + (F4*VECT(IPZ) + F5*HXP(3) + F6*H(3))
 
      GO TO 999
 
   10 CONTINUE
      DO 20 I   = 1,3
         VOUT(I) = VECT(I) + STEP * VECT(I+3)
         VOUT(I+3) = VECT(I+3)
   20 CONTINUE
C
  999 END
