*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.41  by  S.Giani
*-- Author :
*$ CREATE CALUMV.FOR
*COPY CALUMV
*
*=== calumv ===========================================================*
*
      SUBROUTINE CALUMV (N,ITTA)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/reac.inc"
#include "geant321/redver.inc"
#include "geant321/split.inc"
*
C*** C.M.S.-ENERGY AND REACTION CHANNEL THRESHOLD CALCULATION
****  INTEGER * 2 ICH,IBAR,K1,K2,NZK,NRK
**** *,IEII,IKII,NURE
      COMMON / FKABLT / AM   (110), GA   (110), TAU  (110), ICH   (110),
     &                  IBAR (110), K1   (110), K2   (110)
      IRE=NURE(N,ITTA/8+1)
      IEO=IEII (IRE)+1
      IEE=IEII (IRE+1)
      AM1  = AM (N)
      AM12 = AM1**2
      AM2  = AM (ITTA)
      AM22 = AM2**2
*  +-------------------------------------------------------------------*
*  |  Loop on the momentum/invariant mass tabulation
      DO 1 IE=IEO,IEE
         PLAB2=PLABF(IE)**2
         ELAB=SQRT(AM12+AM22+2.D+00*SQRT(PLAB2+AM12)*AM2)
         UMO(IE)=ELAB
    1 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      IKO=IKII(IRE)+1
      IKE=IKII(IRE +1)
      UMOO=UMO(IEO)
*  +-------------------------------------------------------------------*
*  |  Loop on the exit channels
      DO 2 IK=IKO,IKE
         INRK1=NRK(1,IK)
         INRK2=NRK(2,IK)
*  |  +----------------------------------------------------------------*
*  |  |  One resonance exit channel
         IF (INRK2.LE.0) THEN
            IF ( GA (INRK1) .GT. ANGLGB ) THEN
               AM111 = AM(INRK1) - 5.D+00 * GA (INRK1)
            ELSE
               AM111 = 0.D+00
            END IF
            INRKK=K1(INRK1)
            AMSS =5.D+00
            INRKO=K2(INRK1)
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Loop over the decay channels
            DO 3 INKK=INRKK,INRKO
               INZK1=NZK(INKK,1)
               INZK2=NZK(INKK,2)
               INZK3=NZK(INKK,3)
*  |  |  |  We allow for a width of at least 1 gamma for resonances !!!
               AMS=AM(INZK1)+AM(INZK2)-GA(INZK1)-GA(INZK2)
               IF (INZK3 .GT. 0) AMS=AMS+AM(INZK3)-GA(INZK3)
               IF (AMSS  .GT.AMS) AMSS=AMS
    3       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            AMS=MAX(AMSS,AM111)
            IF (AMS.LT.UMO(IEO)) AMS=UMO(IEO)
            THRESH(IK)=AMS
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
    2 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      RETURN
      END
