*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE TWOBOD(D,LD,KZ1,KZ2,A1,A2,Z1,Z2,ATAR,Q,MT)
C       THIS ROUTINE CALCULATES THE EXIT ENERGIES AND DIRECTIONAL
C       COSINES FOR THE CHARGED PARTICLE AND RECOIL NUCLEUS FOR
C       A TWO-BODY REACTION USING AN EVAPORATION SPECTRUM AND
C       MOMEMTUM BALANCE.  IT ALSO SETS ALL EXIT PARAMETERS FOR
C       THE COLLISION PRODUCTS AND STORES THEM IN THE RECOIL BANK.
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mnutrn.inc"
#include "geant321/mrecoi.inc"
#include "geant321/mapoll.inc"
#include "geant321/mmass.inc"
#include "geant321/mpstor.inc"
      DIMENSION D(*),LD(*)
      SAVE
      PRXO = 0.0
      PRYO = 0.0
      PRZO = 0.0
C loop over no. of emmitted particles CZ July 30,1992
      NPN = 1
      IF(MT.EQ.108) NPN = 2
      IF(MT.EQ.109) NPN = 3
      IF(MT.EQ.111) NPN = 2
C       CALCULATE THE COULOMB BARRIER (CB)
      CALL BARIER(KZ1,KZ2,A1,A2,CB)
C       CALCULATE THE CHARGED PARTICLE EXIT ENERGY (EX)
      CALL CEVAP(EOLD,Q,ATAR,CB,EX)
      E1=EX+CB
C calculate the massnumber and mass of the residual nucleus
      A2 = A2 - (NPN-1)*A1
      Z2 = Z2 - (NPN-1)*Z1
      IF(A2.LT.0.) A2 = 0.
      IF(Z2.LT.0.) Z2 = 0.
      IF(NPN.EQ.1) THEN
C for 1 final state particle the available kinetic energy is given
C by momentum and energy conservation
        E1 = E1*Z2/(Z1+Z2)
      ENDIF
      DO 10  NP=1,NPN
C       ASSUME ISOTROPIC CHARGED PARTICLE EMISSION IN THE LABORATORY
         CALL GTISO(U1,V1,W1)
         IF(NPN.EQ.1) THEN
C only one final state particle -> use all the energy available
           PPN = SQRT(2.0*Z1*E1)
           EKN = E1
         ELSE
           IF(NP.LT.NPN) THEN
             EKN = E1*FLTRNF(0)
           ELSE
             EKN = E1
           ENDIF
           E1 = E1 - EKN
           PPN = SQRT(2.0*Z1*EKN)
         ENDIF
         PRXO = PRXO + U1*PPN
         PRYO = PRYO + V1*PPN
         PRZO = PRZO + W1*PPN
C       CALCULATE AND SET THE CHARGED PARTICLE EXIT PARAMETERS
         XR=X
         YR=Y
         ZR=Z
         WATER=WTBC
         NZR=KZ1
         AGER=AGE
         NCOLR=NCOL
         MTNR=MT
         AR=A1
         ENIR=EOLD
         UNIR=UOLD
         VNIR=VOLD
         WNIR=WOLD
         ENOR=0.0
         UNOR=0.0
         VNOR=0.0
         WNOR=0.0
         WTNR=0.0
         QR=Q
         UR=U1
         VR=V1
         WR=W1
         ER=EKN
C       STORE THE CHARGED PARTICLE IN THE RECOIL BANK
         EP = ER
         UP = UR
         VP = VR
         WP = WR
         AMP = AR
         ZMP = FLOAT(NZR)
         AGEP = AGE
         MTP = MT
         CALL STOPAR(IDHEVY,NHEVY)
   10 CONTINUE
C       CALCULATE THE TOTAL MOMENTUM BEFORE THE COLLISION
C       NEUTRON MOMENTUM BEFORE COLLISION (PI) EQUALS TOTAL MOMENTUM
      PI=SQRT(2.0*ZN*EOLD)
C       CALCULATE THE DIRECTIONAL MOMENTUM OF THE RECOIL NUCLEUS
      PRX=PI*UOLD - PRXO
      PRY=PI*VOLD - PRYO
      PRZ=PI*WOLD - PRZO
C       CALCULATE THE TOTAL MOMENTUM OF THE RECOIL NUCLEUS
      PR=SQRT(PRX**2+PRY**2+PRZ**2)
C       CALCULATE THE RECOIL NUCLEUS DIRECTIONAL COSINES
      U2=PRX/PR
      V2=PRY/PR
      W2=PRZ/PR
C       CALCULATE THE RECOIL NUCLEUS EXIT ENERGY
      XM  = A2*931.075E6
      E2 = SQRT(PR**2+XM**2) - XM
C       CALCULATE AND SET THE CHARGED PARTICLE EXIT PARAMETERS
      XR=X
      YR=Y
      ZR=Z
      WATER=WTBC
      NZR=KZ2
      AGER=AGE
      NCOLR=NCOL
      MTNR=MT
      AR=A2
      ENIR=EOLD
      UNIR=UOLD
      VNIR=VOLD
      WNIR=WOLD
      ENOR=0.0
      UNOR=0.0
      VNOR=0.0
      WNOR=0.0
      WTNR=0.0
      QR=Q
      UR=U2
      VR=V2
      WR=W2
      ER=E2
C       STORE THE RECOIL HEAVY ION IN THE RECOIL BANK
      EP = ER
      UP = UR
      VP = VR
      WP = WR
      AMP = AR
      ZMP = FLOAT(NZR)
      AGEP = AGE
      MTP = MT
      CALL STOPAR(IDHEVY,NHEVY)
      RETURN
      END
