*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE RCHANV.FOR
*COPY RCHANV
*
*=== rchanv ===========================================================*
*
      SUBROUTINE RCHANV
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/hadflg.inc"
#include "geant321/reac.inc"
#include "geant321/redver.inc"
#include "geant321/split.inc"
*
      COMMON / FKABLT / AM   (110), GA   (110), TAU  (110), ICH   (110),
     &                  IBAR (110), K1   (110), K2   (110)
*  Note: 296 are the total number of energy at which data are tabulated
*        (of course for the 17 reactions considered, depending on the
*         reaction there could be different numbers of tabulated
*         energies)
*        268 is the number of possible exit channels
*        Umo (ieii(ire)+ie) is the cms energy corresponding to the ieth
*        energy at which data are tabulated for the reaction ire
*        Plabf (ieii(ire)+ie) is the corresponding lab momentum
*        Siin  (ieii(ire)+ie) is the cross section
      DIMENSION HWT(460)
      DIMENSION HWK(40)
      DIMENSION SI(5184)
      EQUIVALENCE (WK(1),SI(1))
C*** WEIGHTS FOR THE SAMPLING PROCEDURE (ADDED ONE TO EACH OTHER IN
C*** CORRESP. CHANNELS) SPECIFIC FOR NUCRIN ONLY
C*** CALCULATION OF THRESHOLD ENERGY OF THE REACTION CHANNELS
C
      IREG=16
*  +-------------------------------------------------------------------*
*  |  Loop on the possible reactions (pi+ p, .... )
      DO 222 IRE=1,IREG
*  |  Initial index for the exit channel sigmas/weights for reaction IRE
*  |  (wk(ire+1)-wk(ire+iee), weights at the various energies for the
*  |  first channel, wk(ire+(ik-1)*iee+ie), weight of the ikth channel
*  |  at ieth energy)
         IWKO=IRII(IRE)
*  |  Number of energy tabulations for reaction ire
         IEE=IEII(IRE+1)-IEII(IRE)
*  |  Number of exit channels of reaction ire
         IKE=IKII(IRE+1)-IKII(IRE)
*  |  Index for the initial energy tabulation for reaction ire (this is
*  |  for index 1!!, ieii is for index 0)
         IEO=IEII(IRE)+1
*  |  Index for the initial exit channel of the reaction ire
*  |  (the initial channel is IIKI + 1)
         IIKI=IKII(IRE)
*  |  +----------------------------------------------------------------*
*  |  |  This loop checks the threshold (expressed in invariant mass)
*  |  |  for the several reaction channels:
*  |  |  Channels resulting in two exit particles/resonances are
*  |  |  checked for Thresh >= m(1) + m(2)
*  |  |  Channels resulting in only one resonance are checked for
*  |  |  Thresh >= Min_j (m_j(1)+m_j(2)+m_j(3)), where the minimum
*  |  |  is carried out looping over all possible decay channels j
*  |  |  and now also looking for the mass of the resonance
*  |  |  less 5 x width
         DO 226 IK=1,IKE
            INRK1 = NRK(1,IIKI+IK)
            INRK2 = NRK(2,IIKI+IK)
            AM111 = AM (INRK1)
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Two particles/resonances exit channels
            IF ( INRK2 .GT. 0 ) THEN
               AM222 = AM (INRK2)
               THRESH (IIKI+IK) = AM111 + AM222
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  One resonance exit channel
            ELSE
               IF ( GA (INRK1) .GT. ANGLGB ) THEN
                  AM111 = AM111 - 5.D+00 * GA (INRK1)
               ELSE
                  AM111 = 0.D+00
               END IF
               INRKK = K1(INRK1)
               AMSS  = 5.D+00
               INRKO = K2(INRK1)
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Loop over the decay channels
               DO 228 INKK=INRKK,INRKO
                  INZK1=NZK(INKK,1)
                  INZK2=NZK(INKK,2)
                  INZK3=NZK(INKK,3)
                  AMS = AM(INZK1)+AM(INZK2)-2.D+00*(GA(INZK1)+GA(INZK2))
                  IF (INZK3 .GT. 0)  AMS =AMS+AM(INZK3)-2.D+00*GA(INZK3)
                  IF (AMSS  .GT.AMS) AMSS=AMS
  228          CONTINUE
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               AMS = MAX (AMSS,AM111)
               IF ( AMS .LT. UMO(IEO) ) AMS = UMO (IEO)
               THRESH (IIKI+IK) = AMS
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
  226    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         SINORC = 1.D+00
*  |  +----------------------------------------------------------------*
*  |  |  Loop on the energy tabulations
         DO 221 IE=1,IEE
            SIS=ANGLGB/10.D+00
            PLASQ = PLABF (IEO+IE-1)**2
            UMOSQ = ( SQRT ( AM (INNURE(1,1,IRE))**2 + PLASQ )
     &            + AM (INNURE(2,1,IRE)) )**2 - PLASQ
            IF ( INNURE (1,2,IRE) .GT. 0 )
     &         UMOSQ = MAX ( UMOSQ, ( SQRT ( AM (INNURE(1,2,IRE))**2
     &               + PLASQ ) + AM (INNURE(2,2,IRE)) )**2 - PLASQ )
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Loop on the exit channels
            DO 223 IK=1,IKE
*  |  |  |  IWK index of the sigma (weight) of the IKth exit channel of
*  |  |  |  reaction IRE at energy IE
               IWK=IWKO+IEE*(IK-1)+IE
*  |  |  |  NRK (i,iiki+ik), i=1,2 are the two resonances produced by
*  |  |  |  the exit channel ik of the reaction ire: 0 means no second
*  |  |  |  resonance
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Check that cross section is 0 below the computed
*  |  |  |  |  threshold
               IF ( UMOSQ .GE. THRESH (IIKI+IK)**2 ) THEN
                  SIS=SIS+SI(IWK)*SINORC
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
                  SI(IWK)=0.D+00
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
  223       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            SIIN(IEO+IE-1)=SIS
            SIO=0.D+00
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF (SIS.LE.ANGLGB) THEN
               SIS=1.D+00
               SIO=1.D+00
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 224 IK=1,IKE
               IWK=IWKO+IEE*(IK-1)+IE
               SIO=SIO+SI(IWK)/SIS*SINORC
               HWK(IK)=SIO
  224       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 225 IK=1,IKE
               IWK=IWKO+IEE*(IK-1)+IE
               WK(IWK)=HWK(IK)
  225       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
  221    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
  222 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      DO 3 J=1,460
         HWT(J)=0.D+00
    3 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      DO 1 I=1,110
         IK1=K1(I)
         IK2=K2(I)
         HV=0.D+00
         DO 2 J=IK1,IK2
            HV=HV+WT(J)
            HWT(J)=HV
            JI=J
    2    CONTINUE
         IF (ABS(HV-1.D0).GT.1.D-4)WRITE(LUNOUT,101)
  101    FORMAT(44H ERROR IN HWT BECAUSE OF FALSE USE OF RCHANW)
    1 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      DO 4 J=1,460
         WT(J)=HWT(J)
    4 CONTINUE
*  |
*  +-------------------------------------------------------------------*
* Set a flag for hadrin that elastic collisions must be reduced
* because they will occur inside nuclei
      IELFLG = -1
      ICXFLG = -1
      RETURN
      END
