*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDECA3(XM0,XM1,XM2,XM3,PCM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   This routine simulates three-body decays in center-of-mass.  *
C.    *   Written by Jurgen Schukraft and Vincent Hedberg for R807.    *
C.    *   Adapted for Geant3 by Sverker Johansson (Aug-1987).          *
C.    *   (Quite extensive modifications; comments etc. are not always *
C.    *    consistent with the code anymore. )                         *
C.    *                                                                *
C.    *   As it is, the decay is done according to pure phase-space.   *
C.    *   There is code in the routine to properly treat K decays (V-A)*
C.    *   but it is not used. (Requires an additional parameter        *
C.    *   in the call, to tell which type of decay.)                   *
C.    *                                                                *
C.    *    ==>Called by : GDECAY                                       *
C.    *       Author    S. Johansson  *********                        *
C.    *                                                                *
C.    ******************************************************************
C.
***************************************************************
******* original name :                     *******************
******SUBROUTINE  D E C 3 B (ID,IP,K1,K2,K3)*******************
***************************************************************
C--------------------------------------------------------
C      O B S O L E T E   C O M M E N T   ! ! !
C     THIS ROUTINE SIMULATES THE THREE-BODY DECAYS OF:
C--------------------------------------------------------
C     ID = 1 : ETA   INTO  PI0 + PI0 + PI0     PHASE SPACE
C     ID = 2 : OMEGA INTO  PI0 + PI- + PI+        - " -
C     ID = 3 : ETAPR INTO  PI0 + PI0 + ETA        - " -
C     ID = 4 : K+/-  INTO  E+/- + PI0 + NU        (V-A)
C     ID = 5 : K0L   INTO  E+/- + PI-/+ + NU      - " -
C     ID = 6 : D     INTO  E   +  K  + NU      PHASE SPACE
C     ID = 7 : D     INTO  E   +  K* + NU         - " -
**      (** or any other decay; input is now masses
**          rather than process ID   S.J.
**          But the ID parameter could someday be useful
**          since it controls the non-phase-space code. ***)
C--------------------------------------------------------
C     THE OUTPUT COMMON BLOCK /LAB/ PLAB(20,J) IS FILLED
*       (** Not a common block anymore.  S.J. **)
C     ACCORDING TO THE ORDER OF THE SECONDARY PARTICLES
C     AS THEY APPEARE IN THE COMMENTS ABOVE. FOR EXAMPLE,
C     IN THE CASE OF ID=6 (D-MESON DECAY)
C     J = K1 IN PLAB(20,J) CORRESPONDS TO ELECTRON,
C     J = K2  -- " -- " -- " -- " --  K-MESON,
C     J = K3  -- " -- " -- " -- " --  NEUTRINO
C--------------------------------------------------------
      DIMENSION  P(20,20)
      DIMENSION  HM(4,7)
      DIMENSION  PCM(4,3)
      DIMENSION RNDM(3)
*
****************************************************
*   Translation from new to old input parameters :
*      ( S.J. )
****************************************************
      ID = 1
      HM(1,ID) = XM1
      HM(2,ID) = XM2
      HM(3,ID) = XM3
      HM(4,ID) = XM0
      K1 = 1
      K2 = 2
      K3 = 3
      PI2 = 2. * 3.141592654
***********************
*  Original code :
C--------------------------------------------------------
C     SIMULATION OF SECONDARY PARTICLE MOMENTA
C     IN THE REST FRAME OF PARENT PARTICLE
C--------------------------------------------------------
      T0=HM(4,ID)-HM(1,ID)-HM(2,ID)-HM(3,ID)
  10  CALL GRNDM(RNDM,2)
      G1=RNDM(1)
      G2=RNDM(2)
      GMI=MIN(G1,G2)
      GMA=MAX(G1,G2)
      TE=GMI*T0
      TM=(1.-GMA)*T0
      TN=(GMA-GMI)*T0
      PA1=SQRT(TE**2+2.*HM(1,ID)*TE)
      PA2=SQRT(TM**2+2.*HM(2,ID)*TM)
      PA3=SQRT(TN**2+2.*HM(3,ID)*TN)
C--------------------------------------------------------
C     CHECK OF THE MOMENTUM CONSERVATION
C--------------------------------------------------------
      PMMM=MAX(PA1,PA2)
      PMAX=MAX(PA3,PMMM)
      SP=PA1+PA2+PA3
      IF(PMAX.GE.SP-PMAX) GOTO 10
      P(4,K1)=TE+HM(1,ID)
      P(4,K2)=TM+HM(2,ID)
      P(4,K3)=TN+HM(3,ID)
C--------------------------------------------------------
C     ACCOUNTING FOR MATRIX ELEMENT (FOR ID=4 AND 5)
*     (** Calculates correct (V-A) decay for kaons.       **)
*     (** Now inactive.  Change ID to activate it. (S.J)  **)
C--------------------------------------------------------
      IF(ID.NE.4.AND.ID.NE.5) GOTO 20
      EEM=((HM(4,ID)-HM(2,ID))**2+HM(1,ID)**2-
     *      HM(3,ID)**2)/(2.*(HM(4,ID)-HM(2,ID)))
      HMMAX=EEM*(HM(4,ID)-HM(2,ID)-EEM)+
     *      0.5*(EEM**2-HM(1,ID)**2)
      HMA=P(4,K1)*P(4,K3)+0.25*(PA1**2+PA3**2-PA2**2)
      CALL GRNDM(RNDM,1)
      GG=RNDM(1)
      IF(HMA.LT.GG*HMMAX) GOTO 10
C--------------------------------------------------------
C     CALCULATIONS OF MOMENTA
C--------------------------------------------------------
20    CONTINUE
      CALL GRNDM(RNDM,3)
      CTE=2.*RNDM(1)-1.
      STE=SQRT(ABS(1.-CTE**2))
      FE =PI2*RNDM(2)
      CFE=COS(FE)
      SFE=SIN(FE)
      P(1,K1)=PA1*STE*CFE
      P(2,K1)=PA1*STE*SFE
      P(3,K1)=PA1*CTE
      CTEN=(PA2**2-PA1**2-PA3**2)/(2.*PA1*PA3)
      STEN=SQRT(ABS(1.-CTEN**2))
      FEN =PI2*RNDM(3)
      CFEN=COS(FEN)
      SFEN=SIN(FEN)
      P(1,K3)=PA3*(STEN*CFEN*CTE*CFE-STEN*SFEN*SFE+CTEN*STE*CFE)
      P(2,K3)=PA3*(STEN*CFEN*CTE*SFE+STEN*SFEN*CFE+CTEN*STE*SFE)
      P(3,K3)=PA3*(-STEN*CFEN*STE+CTEN*CTE)
      P(1,K2)=-P(1,K1)-P(1,K3)
      P(2,K2)=-P(2,K1)-P(2,K3)
      P(3,K2)=-P(3,K1)-P(3,K3)
C-------------------------------------------------------
***   Re-translation from old to new output arrays : ***
***     (S.J.)                                       ***
***--------------------------------------------------***
      DO 37 I=1,4
        DO 38 J=1,3
          PCM(I,J) = P(I,J)
 38     CONTINUE
 37   CONTINUE
      END
