*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:34  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERLAND (STEP, Z, A, RHO, P, E, XMASS, DEDX2)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates energy straggling using Gaussian theory in a step  *
C.    *                                                                *
C.    *  Input  STEP   =  current step-length (cm)                     *
C.    *         Z      =  Z of the material                            *
C.    *         A      =  A of the material                            *
C.    *         RHO    =  density of the material                      *
C.    *         P      =  momentum of the particle                     *
C.    *         E      =  energy   of the particle                     *
C.    *         XMASS  =  mass     of the particle                     *
C.    *                                                                *
C.    *  Output DEDX2  =  mean square of the straggling in G           *
C.    *                                                                *
C.    *    ==>Called by : ERTRCH                                       *
C.    *       Author    E.Nagy  *********                              *
C.    *       Original routine : GLANDO                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconst.inc"
*
      IF (STEP.LT.1.E-7) THEN
         DEDX2=0.
         RETURN
      ENDIF
*
*     Calculate xi factor (KeV).
      BETA   = P/E
      GAMMA  = E/XMASS
      XI     = (153.5*Z*STEP*RHO)/(A*BETA*BETA)
*
*     Maximum energy transfer to atomic electron (KeV).
      ETA    = BETA*GAMMA
      ETASQ  = ETA*ETA
      RATIO  = EMASS/XMASS
      F1     = 2.*EMASS*ETASQ
      F2     = 1.+2.*RATIO*GAMMA+RATIO*RATIO
      EMAX   = 1.E+6*F1/F2
*
*
* *** and now sigma**2  in GeV
*        gaussian
      DEDX2  = XI*EMAX*(1.-(BETA*BETA/2.))*1.E-12
*
      END
