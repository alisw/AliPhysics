*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDRELI(A1,Z1,A2,Z2,DENS,T,DEDX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates the mean 1/DENS*dE/dx of an ion with charge Z1,    *
C.    *  atomic weight A1 and kinetic energy T in an element           *
C.    *  of atomic number Z2, atomic weight A2 and density             *
C.    *    DENS  ( the density is just used for the calculation        *
C.    *  of the density effect in the case of high T ).                *
C.    *  The routine reproduces the experimental and/or tabulated      *
C.    *  mean energy losses reasonably well.                           *
C.    *                                                                *
C.    *  The mean stopping power is obtained by calculating the proton *
C.    *  energy loss S at the equivalent proton kinetic energy  and    *
C.    *  multiplying this value by the effective charge of the ion.    *
C.    *  This method is used for high T ( T/A1 > Tamax , where         *
C.    *   Tamax depend on Z1 and Z2 , Tamax .le. few MeV).             *
C.    *  In the case of low T , the energy loss curve has been         *
C.    *  extrapolated down to T -> 0.                                  *
C.    *                                                                *
C.    *    ==>Called by : GDRELA                                       *
C.    *       Author    L.Urban    *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcunit.inc"
      PARAMETER (AMU=0.93149432,D=0.00015357)
*
      EM=A1*AMU
      TA=PMASS*T/EM
*
      Z11=Z1*Z1
      TAMAX=(120.902-3.121*Z1+0.270*Z11)-(35.988-27.794*Z1+0.120*Z11)*
     +                                                   LOG(Z2)
      TAMAX=0.000001*TAMAX
*
*     calculate  stopping power (total loss) , save DCUTM before
*
      CUTSAV=DCUTM
      DCUTM=BIG
*
      IF(TA.GE.TAMAX) THEN
*
*     "high energy" case
*
         CALL GDRELP(A2,Z2,DENS,TA,S)
         S=Z1**2*GEFCH2(Z1,Z2,TA)*S
      ELSE
*
*     "low energy" case
*
         CALL GDRELP(A2,Z2,DENS,TAMAX,S0)
         S0=Z1**2*GEFCH2(Z1,Z2,TAMAX)*S0
         R=TA/TAMAX
         S=S0*(2.*SQRT(R)-R)
      ENDIF
*
      DCUTM=CUTSAV
      TMAX=2.*EMASS*T*(T+2.*EM)/EM**2
*
*     check whether restricted loss needed ?
*     if restricted loss , calculate the loss from delta rays
*
      IF(DCUTM.LT.TMAX) THEN
         BET2=T*(T+2.*EM)/(T+EM)**2
         R=DCUTM/TMAX
         DELTA=-LOG(R)-BET2*(1.-R)
         DELTA=D*Z2*Z1**2*DELTA/(A2*BET2)
         S=S-DELTA
         IF(S.LT.0.) S=0.
      ENDIF
*
      DEDX=S
*
      END
