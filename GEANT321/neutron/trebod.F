*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE TREBOD(D,LD,KZ1,KZ2,KZ3,A1,A2,A3,Z1,Z2,Z3,
     +                  ATAR,Q,MT)
C CZ July 30,1992 Simple aproach to get (N,PA), (N,T2A),(N,D2A)
C processes. This is TWOBOD extended to a third particle
C       THIS ROUTINE CALCULATES THE EXIT ENERGIES AND DIRECTIONAL
C       COSINES FOR THE CHARGED PARTICLE AND RECOIL NUCLEUS FOR
C       A THREE-BODY REACTION USING AN EVAPORATION SPECTRUM AND
C       MOMEMTUM BALANCE.  IT ALSO SETS ALL EXIT PARAMETERS FOR
C       THE COLLISION PRODUCTS AND STORES THEM IN THE RECOIL BANK.
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mnutrn.inc"
#include "geant321/mrecoi.inc"
#include "geant321/mapoll.inc"
#include "geant321/mmass.inc"
#include "geant321/mpstor.inc"
      DIMENSION D(*),LD(*),ER1(3)
      SAVE
C loop over no. of emmitted particles CZ July 30,1992
      NPN = 1
      IF(MT.EQ.112) NPN = 2
      IF(MT.EQ.113) NPN = 3
      IF(MT.EQ.114) NPN = 3
      PRXO = 0.0
      PRYO = 0.0
      PRZO = 0.0
      DO 10  NP=1,NPN
C       CALCULATE THE COULOMB BARRIER (CB)
         CALL BARIER(KZ1,KZ2,A1,A3,CB)
C       CALCULATE THE CHARGED PARTICLE EXIT ENERGY (EX)
         CALL CEVAP(EOLD,Q,ATAR,CB,EX)
         E1=EX+CB
         ZMSS = Z2
         AMSS = A2
         KZZ = KZ2
         IF(NP.EQ.1) THEN
            ZMSS = Z1
            AMSS = A1
            KZZ = KZ1
         ENDIF
C       ASSUME ISOTROPIC CHARGED PARTICLE EMISSION IN THE LABORATORY
         CALL GTISO(U1,V1,W1)
         PPN = SQRT(2.0*ZMSS*E1)
         PRXO = PRXO + U1*PPN
         PRYO = PRYO + V1*PPN
         PRZO = PRZO + W1*PPN
C       CALCULATE AND SET THE CHARGED PARTICLE EXIT PARAMETERS
         XR=X
         YR=Y
         ZR=Z
         WATER=WTBC
         NZR=KZZ
         AGER=AGE
         NCOLR=NCOL
         MTNR=MT
         AR=AMSS
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
         ER=E1
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
         A3 = A3 - A2
         Z3 = Z3 - Z2
         KZ3 = KZ3 - KZ2
   10 CONTINUE
      A3 = A3 + A2
      Z3 = Z3 + Z2
      KZ3 = KZ3 + KZ2
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
      XM = A2 * 931.075E6
      E2 = SQRT(PR**2+XM**2) - XM
C       CALCULATE AND SET THE CHARGED PARTICLE EXIT PARAMETERS
      XR=X
      YR=Y
      ZR=Z
      WATER=WTBC
      NZR=KZ3
      AGER=AGE
      NCOLR=NCOL
      MTNR=MT
      AR=A3
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
