*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:33  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GPROBI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Initialise material constants used in the computation of      *
C.    *  the probability for various interactions.                     *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Authors    R.Brun, G.Patrick, L.Urban  *********         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
      DIMENSION EK(4),EL1(4),EL2(4)
      DATA EK / 0.66644E-8 , 0.22077E-9 ,-0.32552E-11, 0.18199E-13/
      DATA EL1/-0.29179E-9 , 0.87983E-10,-0.12589E-11, 0.69602E-14/
      DATA EL2/-0.68606E-9 , 0.10078E-9 ,-0.14496E-11, 0.78809E-14/
      DATA ALFA/7.29735E-3/
      DATA REL/0.2817938/
C.
C.    ------------------------------------------------------------------
C.
      IF(Z.LT.1.) GOTO 999
      AEFF=A
      JPROB = LQ(JMA-4)
      IF(JMIXT.GT.0)THEN
         JMI1=LQ(JMIXT-1)
         AEFF=Q(JMI1+1)
      ENDIF
C
C             store constants for PAIR/BREMS routines
C
      X      = (Z*ALFA)**2
      FC     = (( - 0.002 * X + 0.0083) * X - 0.0369) * X + 0.20206
      FC     = X * (FC + 1. / (1. + X))
      C1=Z**0.333333
      C2=LOG(C1)
      C3=LOG(183./C1)-FC
      C4=LOG(1440./(C1*C1))/C3
      Q(JPROB+1)=Z*(Z+C4)*C3/A
      Q(JPROB+2)=C1
      Q(JPROB+3)=C2
      Q(JPROB+4)=FC
C
C             constants for PHOTOEFFECT
C
      Z2   = Z*Z
      EKZ  = Z2*(EK(1) +Z*(EK(2) +Z*(EK(3) +Z*EK(4))))
      EL1Z = Z2*(EL1(1)+Z*(EL1(2)+Z*(EL1(3)+Z*EL1(4))))
      EL2Z = Z2*(EL2(1)+Z*(EL2(2)+Z*(EL2(3)+Z*EL2(4))))
      Q(JPROB+5)=EKZ
      Q(JPROB+6)=EL1Z
      Q(JPROB+7)=EL2Z
C
C             Constants for Hadronic interactions
C
      Q(JPROB+8)= 1000.*AEFF/(AVO*DENS)
C
C             Constants for electron/positron ionisation losses
C             and S5 for one-photon annihilation
C
      IF(JMIXT.LE.0)THEN
         POTI=16.E-9*Z**0.9
         S1=Z/A
         S5=Z**5/A*ALFA**4
      ELSE
         NLMAT=Q(JMA+11)
         NLM=IABS(NLMAT)
         S1=0.
         S2=0.
         S5=0.
         DO 10 J=1,NLM
            AJ=Q(JMIXT+J)
            ZJ=Q(JMIXT+NLM+J)
            WJ=Q(JMIXT+2*NLM+J)
            S1=S1+WJ*ZJ/AJ
            S2=S2+WJ*ZJ*LOG(ZJ)/AJ
            S5=S5+WJ*ZJ**5/AJ*ALFA**4
   10    CONTINUE
         POTI=16.E-9*EXP(0.9*S2/S1)
      ENDIF
      Q(JPROB+9) = POTI
      Q(JPROB+10) = LOG(POTI)
C
      CON1=LOG(POTI/EMASS)
      CON2=DENS*S1
      CON3=1.+2.*LOG(POTI/(28.8E-9*SQRT(CON2)))
C
C             Condensed material ?
C             (at present that means: DENS.GT.0.05 g/cm**3)
C
      IF(DENS.GT.0.05)THEN
         IF(POTI.LT.1.E-7)THEN
            IF(CON3.LT.3.681)THEN
               CON4=0.2
            ELSE
               CON4=0.326*CON3-1.
            ENDIF
            CON5=2.
         ELSE
            IF(CON3.LT.5.215)THEN
               CON4=0.2
            ELSE
               CON4=0.326*CON3-1.5
            ENDIF
            CON5=3.
         ENDIF
      ELSE
C
C             Gas (T=0 C, P= 1 ATM)
C             if T.NE. 0 C and/or P.NE. 1 ATM
C             you have to modify the variable X
C             X=>X+0.5*LOG((273+T C)/(273*P ATM))
C             in the function GDRELE
C             ------------------------
C
         IF(CON3.LE.12.25)THEN
            IP=INT((CON3-10.)/0.5)+1
            IF(IP.LT.0) IP=0
            IF(IP.GT.4) IP=4
            CON4=1.6+0.1*FLOAT(IP)
            CON5=4.
         ELSE
            IF(CON3.LE.13.804)THEN
               CON4=2.
               CON5=5.
            ELSE
               CON4=0.326*CON3-2.5
               CON5=5.
            ENDIF
         ENDIF
      ENDIF
C
      XA=CON3/4.606
      CON6=4.606*(XA-CON4)/(CON5-CON4)**3.
      Q(JPROB+11)=CON1
      Q(JPROB+12)=CON2
      Q(JPROB+13)=-CON3
      Q(JPROB+14)=CON4
      Q(JPROB+15)=CON5
      Q(JPROB+16)=CON6
C
C            constant for delta rays
C            (the same constant is used in the Compton
C              and Annihilation subroutines )
C            and for one-photon annihilation
C
      Q(JPROB+17)=AVO*TWOPI*REL*REL*DENS*S1
      Q(JPROB+18)=AVO*TWOPI*REL*REL*DENS*S5
C
C            Constants for Moliere scattering
C
      IF(JMIXT.LE.0)THEN
         CALL GMOLI(A,Z,1.,1,DENS,Q(JPROB+21),Q(JPROB+25))
      ELSE
         CALL GMOLI(Q(JMIXT+1),Q(JMIXT+NLM+1),Q(JMIXT+2*NLM+1),
     +              NLM,DENS,Q(JPROB+21),Q(JPROB+25))
      ENDIF
C
C                Constants for muon bremsstrahlung
C
      Q(JPROB+31)=LOG(189.*EMMU/(EMASS*C1))
      IF(Z.GT.10)Q(JPROB+31)=Q(JPROB+31)+LOG(0.666666/C1)
      SE         =SQRT(2.71828)
      Q(JPROB+32)=189.*SE*EMMU*EMMU/(2.*EMASS*C1)
      Q(JPROB+33)=0.75*SE*EMMU*C1
C
  999 END
