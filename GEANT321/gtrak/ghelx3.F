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
      SUBROUTINE GHELX3 (FIELD, STEP, VECT, VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Tracking routine in a constant field oriented            *
C.    *       along axis 3                                             *
C.    *       Tracking is performed with a conventional                *
C.    *       helix step method                                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Authors    R.Brun, M.Hansroul  *********                 *
C     *       Rewritten  V.Perevoztchikov
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION      VECT(7),VOUT(7),HXP(3)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6,IPP=7)
      PARAMETER      (SIXTH = 1./6.)
      PARAMETER      (EC=2.9979251E-4)
C.
C.    ------------------------------------------------------------------
C.
C       units are kgauss,centimeters,gev/c
C
      VOUT(IPP) = VECT(IPP)
      H4 = FIELD * EC
*
      HXP(1) = - VECT(IPY)
      HXP(2) = + VECT(IPX)
 
      HP = VECT(IPZ)
*
      RHO = -H4/VECT(IPP)
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
 
      VOUT(IX) = VECT(IX) + (F1*VECT(IPX) + F2*HXP(1))
      VOUT(IY) = VECT(IY) + (F1*VECT(IPY) + F2*HXP(2))
      VOUT(IZ) = VECT(IZ) + (F1*VECT(IPZ) + F3)
 
      VOUT(IPX) = VECT(IPX) + (F4*VECT(IPX) + F5*HXP(1))
      VOUT(IPY) = VECT(IPY) + (F4*VECT(IPY) + F5*HXP(2))
      VOUT(IPZ) = VECT(IPZ) + (F4*VECT(IPZ) + F6)
 
C
  999 END
