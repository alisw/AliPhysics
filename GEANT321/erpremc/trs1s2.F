*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:36  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE TRS1S2 (PD1,RD1,PD2,RD2,H,CH,IERR,SP1,SP2
     1,                  DJ1,DK1,DJ2,DK2)
C
C *** TRANSFORMS ERROR MATRIX
C     FROM        VARIABLES (1/P,V1',W1',V1,W1)
C      TO         VARIABLES (1/P,V2',W2',V2,W2)
C
C     Authors: A. Haas and W. Wittek
C
C
C *** PD1(3)    1/P,V1',W1'                             INPUT
C     PD2(3)    1/P,V2',W2'                            OUTPUT
C     H(3)      MAGNETIC FIELD                          INPUT
C     RD1(15)   ERROR MATRIX IN 1/P,V1',W1',V1,W1       INPUT      (TRIANGLE)
C     RD2(15)   ERROR MATRIX IN 1/P,V2',W2',V2,W2      OUTPUT      (TRIANGLE)
C     CH        CHARGE OF PARTICLE                      INPUT
C               CHARGE AND MAGNETIC FIELD ARE NEEDED
C               FOR CORRELATION TERMS (V2',V1),(V2',W1),(W2',V1),(W2',W1)
C               THESE CORRELATION TERMS APPEAR BECAUSE RD1 IS ASSUMED
C               TO BE THE ERROR MATRIX FOR FIXED U1
C               AND RD2 FOR FIXED U2
C     SP1       SIGN OF U1-COMPONENT OF PARTICLE MOMENTUM     INPUT
C     SP2       SIGN OF U2-COMPONENT OF PARTICLE MOMENTUM    OUTPUT
C     DJ1(3)    UNIT VECTOR IN V1-DIRECTION
C     DK1(3)    UNIT VECTOR IN W1-DIRECTION    OF SYSTEM 1
C     DJ2(3)    UNIT VECTOR IN V2-DIRECTION
C     DK2(3)    UNIT VECTOR IN W2-DIRECTION    OF SYSTEM 2
C
C     IERR      = 0    TRANSFORMATION OK
C               = 1    MOMENTUM PERPENDICULAR TO U2-DIRECTION (V2',W2' NOT DEFIN
C               = 2    MOMENTUM PERPENDICULAR TO X-AXIS
C
C
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL   PD1,PD2,RD1,RD2,H,CH,SP1,SP2,DJ1,DK1,DJ2,DK2
#endif
#include "geant321/trcom3.inc"
      DIMENSION PD1(3),PD2(3),RD1(15),RD2(15),H(3),DJ1(3),DK1(3)
     +,DJ2(3),DK2(3),UN(3),VN(3),DI1(3),DI2(3),TVW1(3),TVW2(3)
C
#if defined(CERNLIB_SINGLE)
      DATA CFACT8 / 2.997925 E-4 /
#endif
#if !defined(CERNLIB_SINGLE)
      DATA CFACT8 / 2.997925 D-4 /
#endif
C
      IERR=0
      PM=PD1(1)
      TVW1(1)=1./SQRT(1.+PD1(2)**2+PD1(3)**2)
      IF(SP1.LT.0.) TVW1(1)=-TVW1(1)
      TVW1(2)=PD1(2)*TVW1(1)
      TVW1(3)=PD1(3)*TVW1(1)
C
      DI1(1)=DJ1(2)*DK1(3)-DJ1(3)*DK1(2)
      DI1(2)=DJ1(3)*DK1(1)-DJ1(1)*DK1(3)
      DI1(3)=DJ1(1)*DK1(2)-DJ1(2)*DK1(1)
C
      DO 5 I=1,3
         TN(I)=TVW1(1)*DI1(I)+TVW1(2)*DJ1(I)+TVW1(3)*DK1(I)
    5 CONTINUE
C
      DI2(1)=DJ2(2)*DK2(3)-DJ2(3)*DK2(2)
      DI2(2)=DJ2(3)*DK2(1)-DJ2(1)*DK2(3)
      DI2(3)=DJ2(1)*DK2(2)-DJ2(2)*DK2(1)
C
      TVW2(1)=TN(1)*DI2(1)+TN(2)*DI2(2)+TN(3)*DI2(3)
      TVW2(2)=TN(1)*DJ2(1)+TN(2)*DJ2(2)+TN(3)*DJ2(3)
      TVW2(3)=TN(1)*DK2(1)+TN(2)*DK2(2)+TN(3)*DK2(3)
C
      IF(TVW2(1).EQ.0.) GO TO 901
      TR=1./TVW2(1)
      SP2=1.
      IF(TVW2(1).LT.0.) SP2=-1.
      PD2(1)=PD1(1)
      PD2(2)=TVW2(2)*TR
      PD2(3)=TVW2(3)*TR
C
      COSL=SQRT(ABS(1.-TN(3)**2))
      IF(COSL.EQ.0.) GO TO 902
      COSL1=1./COSL
      UN(1)=-TN(2)*COSL1
      UN(2)=TN(1)*COSL1
      UN(3)=0.
C
      VN(1)=-TN(3)*UN(2)
      VN(2)=TN(3)*UN(1)
      VN(3)=COSL
C
      UJ1=UN(1)*DJ1(1)+UN(2)*DJ1(2)+UN(3)*DJ1(3)
      UK1=UN(1)*DK1(1)+UN(2)*DK1(2)+UN(3)*DK1(3)
      VJ1=VN(1)*DJ1(1)+VN(2)*DJ1(2)+VN(3)*DJ1(3)
      VK1=VN(1)*DK1(1)+VN(2)*DK1(2)+VN(3)*DK1(3)
C
      UJ2=UN(1)*DJ2(1)+UN(2)*DJ2(2)+UN(3)*DJ2(3)
      UK2=UN(1)*DK2(1)+UN(2)*DK2(2)+UN(3)*DK2(3)
      VJ2=VN(1)*DJ2(1)+VN(2)*DJ2(2)+VN(3)*DJ2(3)
      VK2=VN(1)*DK2(1)+VN(2)*DK2(2)+VN(3)*DK2(3)
C
      J=0
      DO 10 I=1,5
         DO 4 K=I,5
            J=J+1
            A(I,K)=0.
            A(K,I)=0.
            S(J)=RD1(J)
    4    CONTINUE
   10 CONTINUE
C
      IF(CH.EQ.0.) GO TO 6
      HA=SQRT(H(1)**2+H(2)**2+H(3)**2)
      HAM=HA*PM
      IF(HAM.EQ.0.) GO TO 6
      HM=CH/HA
C
      Q=-HAM*CFACT8
      TT=-Q*TR**3
      SJ1I2=DJ1(1)*DI2(1)+DJ1(2)*DI2(2)+DJ1(3)*DI2(3)
      SK1I2=DK1(1)*DI2(1)+DK1(2)*DI2(2)+DK1(3)*DI2(3)
      SK2U=DK2(1)*UN(1)+DK2(2)*UN(2)+DK2(3)*UN(3)
      SK2V=DK2(1)*VN(1)+DK2(2)*VN(2)+DK2(3)*VN(3)
      SJ2U=DJ2(1)*UN(1)+DJ2(2)*UN(2)+DJ2(3)*UN(3)
      SJ2V=DJ2(1)*VN(1)+DJ2(2)*VN(2)+DJ2(3)*VN(3)
C
      SINZ=-(H(1)*UN(1)+H(2)*UN(2)+H(3)*UN(3))*HM
      COSZ= (H(1)*VN(1)+H(2)*VN(2)+H(3)*VN(3))*HM
      A(2,4)=-TT*SJ1I2*(SK2U*SINZ-SK2V*COSZ)
      A(2,5)=-TT*SK1I2*(SK2U*SINZ-SK2V*COSZ)
      A(3,4)= TT*SJ1I2*(SJ2U*SINZ-SJ2V*COSZ)
      A(3,5)= TT*SK1I2*(SJ2U*SINZ-SJ2V*COSZ)
C
    6 A(1,1)=1.
      A(4,4)=TR*(UJ1*VK2-VJ1*UK2)
      A(4,5)=TR*(UK1*VK2-VK1*UK2)
      A(5,4)=TR*(VJ1*UJ2-UJ1*VJ2)
      A(5,5)=TR*(VK1*UJ2-UK1*VJ2)
C
      TS=TR*TVW1(1)
      A(2,2)=A(4,4)*TS
      A(2,3)=A(4,5)*TS
      A(3,2)=A(5,4)*TS
      A(3,3)=A(5,5)*TS
C
      CALL SSMT5T(A,S,S)
C
      J=0
      DO 25 I=1,5
         DO 20 K=I,5
            J=J+1
            RD2(J)=S(J)
   20    CONTINUE
   25 CONTINUE
C
      RETURN
C
C *** ERROR EXITS
C
  901 IERR=1
      GO TO 910
  902 IERR=2
  910 RETURN
      END
*
