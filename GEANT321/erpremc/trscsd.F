*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
C
      SUBROUTINE TRSCSD(PC,RC,PD,RD,H,CH,IERR,SPU,DJ,DK)
C
C *** TRANSFORMS ERROR MATRIX
C     FROM   SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
C      TO         VARIABLES (1/P,V',W',V,W)
C
C     Authors: A. Haas and W. Wittek
C
C
C *** PC(3)     1/P,LAMBDA,PHI                          INPUT
C     PD(3)     1/P,V',W'                              OUTPUT
C     H(3)      MAGNETIC FIELD                          INPUT
C     RC(15)    ERROR MATRIX IN   SC   VARIABLES        INPUT     (TRIANGLE)
C     RD(15)    ERROR MATRIX IN 1/P,V',W',V,W          OUTPUT     (TRIANGLE)
C     CH        CHARGE OF PARTICLE                      INPUT
C               CHARGE AND MAGNETIC FIELD ARE NEEDED
C               FOR CORRELATION TERMS (V',YT),(V',ZT),(W',YT),(W',ZT)
C               THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
C               TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
C               AND RD FOR FIXED U
C     DJ(3)     UNIT VECTOR IN V-DIRECTION
C     DK(3)     UNIT VECTOR IN W-DIRECTION    OF DETECTOR SYSTEM
C
C     IERR  =   1       PARTICLE MOVES PERPENDICULAR TO U-AXIS
C                      ( V',W' ARE NOT DEFINED )
C     SPU       SIGN OF U-COMPONENT OF PARTICLE MOMENTUM   OUTPUT
C
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL  PC,PD,H,RC,RD,CH,DJ,DK,SPU
#endif
#include "geant321/trcom3.inc"
      DIMENSION PC(3),PD(3),H(3),RC(15),RD(15),DJ(3),DK(3)
      DIMENSION UN(3),VN(3),DI(3),TVW(3)
C
#if defined(CERNLIB_SINGLE)
      DATA CFACT8 / 2.997925 E-4 /
#endif
#if !defined(CERNLIB_SINGLE)
      DATA CFACT8 / 2.997925 D-4 /
#endif
C
      IERR=0
      PM=PC(1)
      COSL=COS(PC(2))
      SINP=SIN(PC(3))
      COSP=COS(PC(3))
C
      TN(1)=COSL*COSP
      TN(2)=COSL*SINP
      TN(3)=SIN(PC(2))
C
      DI(1)=DJ(2)*DK(3)-DJ(3)*DK(2)
      DI(2)=DJ(3)*DK(1)-DJ(1)*DK(3)
      DI(3)=DJ(1)*DK(2)-DJ(2)*DK(1)
C
      TVW(1)=TN(1)*DI(1)+TN(2)*DI(2)+TN(3)*DI(3)
      SPU=1.
      IF(TVW(1).LT.0.) SPU=-1.
      TVW(2)=TN(1)*DJ(1)+TN(2)*DJ(2)+TN(3)*DJ(3)
      TVW(3)=TN(1)*DK(1)+TN(2)*DK(2)+TN(3)*DK(3)
      IF(TVW(1).EQ.0.) GO TO 901
C
      T1R=1./TVW(1)
      PD(1)=PC(1)
      PD(2)=TVW(2)*T1R
      PD(3)=TVW(3)*T1R
C
      UN(1)=-SINP
      UN(2)=COSP
      UN(3)=0.
C
      VN(1)=-TN(3)*UN(2)
      VN(2)=TN(3)*UN(1)
      VN(3)=COSL
C
      UJ=UN(1)*DJ(1)+UN(2)*DJ(2)+UN(3)*DJ(3)
      UK=UN(1)*DK(1)+UN(2)*DK(2)+UN(3)*DK(3)
      VJ=VN(1)*DJ(1)+VN(2)*DJ(2)+VN(3)*DJ(3)
      VK=VN(1)*DK(1)+VN(2)*DK(2)+VN(3)*DK(3)
C
C
      J=0
 
      DO 10 I=1,5
         DO 5 K=I,5
            J=J+1
            A(I,K)=0.
            A(K,I)=0.
            S(J)=RC(J)
    5    CONTINUE
   10 CONTINUE
C
      IF(CH.EQ.0.) GO TO 6
      HA=SQRT(H(1)**2+H(2)**2+H(3)**2)
      HAM=HA*PM
      IF(HAM.EQ.0.) GO TO 6
      HM=CH/HA
      Q=-HAM*CFACT8
C
C
      SINZ=-(H(1)*UN(1)+H(2)*UN(2)+H(3)*UN(3))*HM
      COSZ= (H(1)*VN(1)+H(2)*VN(2)+H(3)*VN(3))*HM
      T3R=Q*T1R**3
      UI=UN(1)*DI(1)+UN(2)*DI(2)+UN(3)*DI(3)
      VI=VN(1)*DI(1)+VN(2)*DI(2)+VN(3)*DI(3)
      A(2,4)=-UI*(VK*COSZ-UK*SINZ)*T3R
      A(2,5)=-VI*(VK*COSZ-UK*SINZ)*T3R
      A(3,4)= UI*(VJ*COSZ-UJ*SINZ)*T3R
      A(3,5)= VI*(VJ*COSZ-UJ*SINZ)*T3R
C
    6 T2R=T1R**2
C
      A(1,1)=1.
      A(2,2)=-UK*T2R
      A(2,3)=VK*COSL*T2R
      A(3,2)=UJ*T2R
      A(3,3)=-VJ*COSL*T2R
      A(4,4)=VK*T1R
      A(4,5)=-UK*T1R
      A(5,4)=-VJ*T1R
      A(5,5)=UJ*T1R
C
      CALL SSMT5T(A,S,S)
C
      J=0
      DO 25 I=1,5
         DO 20 K=I,5
            J=J+1
            RD(J)=S(J)
   20    CONTINUE
   25 CONTINUE
C
      RETURN
C
C *** ERROR EXITS
C
  901 IERR=1
  910 CONTINUE
C
      RETURN
      END
*
