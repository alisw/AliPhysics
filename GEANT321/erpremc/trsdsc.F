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
C
      SUBROUTINE TRSDSC(PD,RD,PC,RC,H,CH,IERR,SPU,DJ,DK)
C
C *** TRANSFORMS ERROR MATRIX
C     FROM        VARIABLES (1/P,V',W',V,W)
C      TO    SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
C
C     Authors: A. Haas and W. Wittek
C
C
C *** PD(3)     1/P,V',W'                               INPUT
C     PC(3)     1/P,LAMBDA,PHI                         OUTPUT
C     H(3)      MAGNETIC FIELD                          INPUT
C     RD(15)    ERROR MATRIX IN 1/P,V',W',V,W           INPUT      (TRIANGLE)
C     RC(15)    ERROR MATRIX IN   SC   VARIABLES       OUTPUT      (TRIANGLE)
C     CH        CHARGE OF PARTICLE                      INPUT
C               CHARGE AND MAGNETIC FIELD ARE NEEDED
C               FOR CORRELATION TERMS (LAMBDA,V),(LAMBDA,W),(PHI,V),(PHI,W)
C               THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
C               TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
C               AND RD FOR FIXED U
C     DJ(3)     UNIT VECTOR IN V-DIRECTION
C     DK(3)     UNIT VECTOR IN W-DIRECTION    OF DETECTOR SYSTEM
C
C     IERR              NOT USED
C     SPU       SIGN OF U-COMPONENT OF PARTICLE MOMENTUM    INPUT
C
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL   PD,PC,H,RC,RD,CH,DJ,DK,SPU
#endif
#include "geant321/trcom3.inc"
      DIMENSION PD(3),PC(3),H(3),RC(15),RD(15),DJ(3),DK(3)
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
      PM=PD(1)
      TVW(1)=1./SQRT(1.+PD(2)**2+PD(3)**2)
      IF(SPU.LT.0.) TVW(1)=-TVW(1)
      TVW(2)=PD(2)*TVW(1)
      TVW(3)=PD(3)*TVW(1)
C
      DI(1)=DJ(2)*DK(3)-DJ(3)*DK(2)
      DI(2)=DJ(3)*DK(1)-DJ(1)*DK(3)
      DI(3)=DJ(1)*DK(2)-DJ(2)*DK(1)
C
      DO 5 I=1,3
         TN(I)=TVW(1)*DI(I)+TVW(2)*DJ(I)+TVW(3)*DK(I)
    5 CONTINUE
C
      PC(1)=PD(1)
      PC(2)=ASIN(TN(3))
      IF (ABS (TN(1)) .LT. 1.E-30) TN(1) = 1.E-30
      PC(3) = ATAN2 (TN(2),TN(1))
C
      COSL=SQRT(ABS(1.-TN(3)**2))
      IF (COSL .LT. 1.E-30) COSL = 1.E-30
      COSL1=1./COSL
      UN(1)=-TN(2)*COSL1
      UN(2)=TN(1)*COSL1
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
      J=0
      DO 10 I=1,5
         DO 4 K=I,5
            J=J+1
            A(I,K)=0.
            A(K,I)=0.
            S(J)=RD(J)
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
C
      SINZ=-(H(1)*UN(1)+H(2)*UN(2)+H(3)*UN(3))*HM
      COSZ= (H(1)*VN(1)+H(2)*VN(2)+H(3)*VN(3))*HM
      A(2,4)=-Q*TVW(2)*SINZ
      A(2,5)=-Q*TVW(3)*SINZ
      A(3,4)=-Q*TVW(2)*COSZ*COSL1
      A(3,5)=-Q*TVW(3)*COSZ*COSL1
C
    6 A(1,1)=1.
      A(2,2)=TVW(1)*VJ
      A(2,3)=TVW(1)*VK
      A(3,2)=TVW(1)*UJ*COSL1
      A(3,3)=TVW(1)*UK*COSL1
      A(4,4)=UJ
      A(4,5)=UK
      A(5,4)=VJ
      A(5,5)=VK
C
      CALL SSMT5T(A,S,S)
C
      J=0
      DO 25 I=1,5
         DO 20 K=I,5
            J=J+1
            RC(J)=S(J)
   20    CONTINUE
   25 CONTINUE
C
      RETURN
C
C *** ERROR EXITS
C
      END
*
