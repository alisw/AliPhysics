*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE CALUMO.FOR
*COPY CALUMO
*
*=== calumo ===========================================================*
*
      SUBROUTINE CALUMO (N,ITTA)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** C.M.S.-ENERGY AND REACTION CHANNEL THRESHOLD CALCULATION
*
#include "geant321/reac.inc"
#include "geant321/redver.inc"
#include "geant321/split.inc"
*
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
      DO 1 IE=IEO,IEE
         PLAB2=PLABF(IE)**2
         ELAB=SQRT(AM12+AM22+2.D+00*SQRT(PLAB2+AM12)*AM2)
         UMO(IE)=ELAB
    1 CONTINUE
      IKO=IKII(IRE)+1
      IKE=IKII(IRE +1)
      UMOO=UMO(IEO)
      DO 2 IK=IKO,IKE
         IF (NRK(2,IK).GT.0) GO TO 2
         IKI=NRK(1,IK)
         AMSS=5.D+00
         K11=K1(IKI)
         K22=K2(IKI)
         DO 3 IK1=K11,K22
            IN=NZK(IK1,1)
            AMS=AM(IN)
            IN=NZK(IK1,2)
            IF (IN.GT.0) AMS=AMS+AM(IN)
            IN=NZK(IK1,3)
            IF (IN.GT.0) AMS=AMS+AM(IN)
            IF (AMS.LT.AMSS) AMSS=AMS
    3    CONTINUE
         IF ( UMOO .LT. AMSS ) UMOO=AMSS
         THRESH (IK) = UMOO
    2 CONTINUE
      RETURN
      END
