*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE LABNF(D,LD,FM,AWR,KZ,Q,LIFLAG)
C       THIS ROUTINE CALCULATES THE DIRECTIONAL COSINES FOR THE
C       NEUTRON BORN FROM THE FISSION REACTION.  THIS VERSION OF
C       THE PROGRAM WILL TREAT A FISSION REACTION AS A SCATTERING
C       EVENT WITH THE NEUTRON EMERGING WITH A MODIFIED WEIGHT OF
C       WATE*NU-BAR.  NO PROVISIONS ARE MADE AT THIS TIME TO
C       CALCULATE THE FISSION FRAGMENTS PARAMETERS, HOWEVER A HEAVY
C       RECOIL ION WILL BE STORED (FOR ANALYSIS PURPOSES) WITH
C       ENERGY AND DIRECTION COSINES EQUAL TO ZERO.
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mnutrn.inc"
#include "geant321/mrecoi.inc"
#include "geant321/mapoll.inc"
#include "geant321/mmass.inc"
#include "geant321/mpstor.inc"
      DIMENSION D(*),LD(*),FM(*)
      SAVE
      MT=18
C       CALCULATE THE NEUTRON EXIT DIRECTIONAL COSINES
      POX = 0.0
      POY = 0.0
      POZ = 0.0
      DO 40 KN=1,INEU
         IF(LIFLAG.EQ.1) THEN
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
   30       UP=UOLD*FM(KN)-COSETA*SINPSI*STHETA
            VP=VOLD*FM(KN)+UOLD*COSPHI*COSETA*SINPSI-SINPHI* SINPSI*
     +      SINETA
            WP=WOLD*FM(KN)+UOLD*SINPHI*COSETA*SINPSI+COSPHI* SINPSI*
     +      SINETA
            S=1.0/SQRT(UP**2+VP**2+WP**2)
            UP=UP*S
            VP=VP*S
            WP=WP*S
         ENDIF
         AGEP = AGE
         EP = ENE(KN)
C use only first neutron for recoil calculation in order to ensure
C correct recoil nucleus energy spectrum
         IF(KN.EQ.1) THEN
            PP = SQRT(EP**2 + 2.0*EP*ZN)
            POX = POX + PP*UP
            POY = POY + PP*VP
            POZ = POZ + PP*WP
         ENDIF
         MTP = MT
         CALL STOPAR(IDNEU,NNEU)
   40 CONTINUE
C       SET THE HEAVY RECOIL ION PARAMETERS FOR ANALYSIS TAPE
   50 PI=SQRT(2.0*ZN*EOLD)
      PIX=PI*UOLD
      PIY=PI*VOLD
      PIZ=PI*WOLD
      PRX=PIX-POX
      PRY=PIY-POY
      PRZ=PIZ-POZ
      PR=SQRT(PRX**2+PRY**2+PRZ**2)
      UR=PRX/PR
      VR=PRY/PR
      WR=PRZ/PR
      AR=AWR*AN+AN-INEU*AN
      XM=AR*931.075E6
      ER=SQRT(PR**2+XM**2)-XM
      EP = ER
      UP = UR
      VP = VR
      WP = WR
      AGEP = AGE
      MTP = MT
      XR=X
      YR=Y
      ZR=Z
      WATER=WTBC
      NZR=KZ
      AGER=AGE
      NCOLR=NCOL
      MTNR=MT
      AMP = AR
      ZMP = FLOAT(KZ)
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
C       STORE THE RECOIL HEAVY ION IN THE RECOIL BANK
      CALL STOPAR(IDHEVY,NHEVY)
      RETURN
      END
