*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE TCHOIC.FOR
*COPY TCHOIC
*
*=== tchoic ===========================================================*
*
      SUBROUTINE TCHOIC(T,P,PP,E,EE,I,II,N,AM1,AM2)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C     ****************************
C     TCHOIC CALCULATES A RANDOM VALUE
C     FOR THE FOUR-MOMENTUM-TRANSFER T
C     ****************************
      COMMON /FKGAMR/ REDU,AMO,AMM(15 )
      COMMON/FKABLT/AM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON/FKSLOP/SM(25),BBM(25),BBB(25)
      REAL RNDM(1)
      AMA=AM1
      AMB=AM2
      IF (I.GT.30.AND.II.GT.30) GO TO 1
      III=II
      AM3=AM2
      IF (I.LE.30) GO TO 11
      III=I
      AM3=AM1
   11 CONTINUE
      GO TO 10
    1 CONTINUE
      III=II
      AM3=AM2
      IF (AMA.LE.AMB) GO TO 10
      III=I
      AM3=AM1
   10 CONTINUE
      IB=IBAR(III)
      AMA=AM3
      K=(AMA-0.75D0)/0.05D0
      IF (K-2.LT.0) K=1
      IF (K-26.GE.0) K=25
      IF (IB)21,22,21
   22 BM=BBM(K)
      GO TO 23
   21 BM=BBB(K)
   23 CONTINUE
C     NORMALIZATION
      TMIN=-2.D0*(E*EE-P*PP)+AM(N)**2+AM1  **2
      TMAX=-2.D0*(E*EE+P*PP)+AM(N)**2+AM1  **2
      CALL GRNDM(RNDM,1)
      VB = RNDM(1)
      IF (VB.LT.0.2D0) BM=BM*0.1D0
     **0.5D0
      TMI=BM*TMIN
      TMA=BM*TMAX
      ETMA=0.D0
      IF (ABS(TMA).GT.120.D0) GO TO 100
      ETMA=EXP(TMA)
  100 CONTINUE
      AN=(1.D0/BM)*(EXP(TMI)-ETMA)
C*** RANDOM CHOICE OF THE T - VALUE
      CALL GRNDM(RNDM,1)
      R = RNDM(1)
      T=(1.D0/BM)*LOG(ETMA+R*AN*BM)
      RETURN
      END
