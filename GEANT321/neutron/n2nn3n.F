*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE N2NN3N(D,LD,AWR,KZ,ID,FM,Q,IFLG)
C       THIS ROUTINE CALCULATES THE DIRECTIONAL COSINES FOR THE
C       NEUTRON AND RECOIL NUCLEUS FOR AN N2N OR N3N REACTION
C       USING THE ONE NEUTRON EMMISION MODEL.  IT ALSO SETS ALL
C       EXIT PARAMETRS FOR THE RECOIL NUCLEUS.
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mnutrn.inc"
#include "geant321/mrecoi.inc"
#include "geant321/mapoll.inc"
#include "geant321/mmass.inc"
#include "geant321/mpstor.inc"
      DIMENSION D(*),LD(*),FM(*)
      SAVE
      MT=0
      IF(ID.EQ.8)MT=16
      IF(ID.EQ.9)MT=17
      IF(ID.EQ.12)MT=24
C       IFLG EQUAL TO ONE IMPLIES THE DIRECTION COSINES WERE
C       SELECTED ISOTROPICALLY IN THE LABORATORY COORDINATE SYSTEM
C       CALCULATE THE NEUTRON EXIT DIRECTIONAL COSINES
      POX = 0.0
      POY = 0.0
      POZ = 0.0
      DO 40 KN=1,INEU
         IF(IFLG.EQ.1) THEN
            CALL GTISO(UP,VP,WP)
         ELSE
            SINPSI=SQRT(1.0-FM(KN)**2)
            CALL AZIRN(SINETA,COSETA)
            STHETA=1.0-UOLD**2
            IF(STHETA)20,20,10
   10       STHETA=SQRT(STHETA)
            COSPHI=VOLD/STHETA
            SINPHI=WOLD/STHETA
            GO TO 30
   20       COSPHI=1.0
            SINPHI=0.0
            STHETA=0.0
   30       UP = UOLD*FM(KN)-COSETA*SINPSI*STHETA
            VP = VOLD*FM(KN)+UOLD*COSPHI*COSETA*SINPSI-SINPHI* SINPSI*
     +      SINETA
            WP = WOLD*FM(KN)+UOLD*SINPHI*COSETA*SINPSI+COSPHI* SINPSI*
     +      SINETA
            S=1.0/SQRT(UP*UP+VP*VP+WP*WP)
            UP=UP*S
            VP=VP*S
            WP=WP*S
         ENDIF
         EP = ENE(KN)
C use ONLY first neutron for recoil calculation in order the ensure
C correct energy spectrum of recoil nucleus
         IF(KN.EQ.1) THEN
            PP = SQRT(EP**2 + 2.0*EP*ZN)
            POX = POX + PP*UP
            POY = POY + PP*VP
            POZ = POZ + PP*WP
         ENDIF
         AGEP = AGE
         MTP = MT
         CALL STOPAR(IDNEU,NNEU)
   40 CONTINUE
C       CALCULATE AND SET THE RECOIL NUCLEUS EXIT PARAMETERS
   50 XR=X
      YR=Y
      ZR=Z
      WATER=WTBC
      NZR=KZ
      ZMP = FLOAT(KZ)
      AGER=AGE
      AGEP = AGE
      NCOLR=NCOL
      MTNR=MT
      MTP = MT
      AR = (AWR*AN) - FLOAT(INEU-1)*AN
      AMP = AR
      ENIR=EOLD
      UNIR=UOLD
      VNIR=VOLD
      WNIR=WOLD
      ENOR=E
      UNOR=U
      VNOR=V
      WNOR=W
      WTNR=WATE
      QR=Q
C       CALCULATE THE NEUTRON MOMENTUM BEFORE AND AFTER COLLISION
C       NEUTRON MOMENTUM BEFORE COLLISION (PI) EQUALS TOTAL MOMENTUM
      PI=SQRT(2.0*ZN*EOLD)
C   CALCULATE THE DIRECTIONAL MOMENTUM OF THE RECOIL NUCLEUS
      PIX=PI*UOLD
      PIY=PI*VOLD
      PIZ=PI*WOLD
      PRX = PIX - POX
      PRY = PIY - POY
      PRZ = PIZ - POZ
C       CALCULATE THE TOTAL MOMENTUM OF THE RECOIL NUCLEUS
      PR=SQRT(PRX**2+PRY**2+PRZ**2)
C       CALCULATE THE RECOIL NUCLEUS DIRECTIONAL COSINES
      UR=PRX/PR
      VR=PRY/PR
      WR=PRZ/PR
      UP = UR
      VP = VR
      WP = WR
C       CALCULATE THE RECOIL NUCLEUS EXIT ENERGY
      XM = AR*931.075E6
      ER= SQRT(PR**2 + XM**2) - XM
      EP = ER
      MTP = MT
C       IF MT=24, DO NOT STORE THE RECOIL HEAVY ION IN THE BANK
      IF(MT.EQ.24)RETURN
C       STORE THE  RECOIL HEAVY ION IN THE RECOIL BANK
      CALL STOPAR(IDHEVY,NHEVY)
      RETURN
      END
