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
      SUBROUTINE TRSPSC(PS,RS,PC,RC,H,CH,IERR,SPX)
C
C *** TRANSFORMS ERROR MATRIX
C     FROM SPLINE VARIABLES (1/P,Y',Z',Y,Z)
C      TO    SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
C
C     Authors: A. Haas and W. Wittek
C
C
C *** PS(3)     1/P,Y',Z'                               INPUT
C     PC(3)     1/P,LAMBDA,PHI                         OUTPUT
C     H(3)      MAGNETIC FIELD                          INPUT
C     RS(15)    ERROR MATRIX IN SPLINE VARIABLES        INPUT      (TRIANGLE)
C     RC(15)    ERROR MATRIX IN   SC   VARIABLES       OUTPUT      (TRIANGLE)
C     CH        CHARGE OF PARTICLE                      INPUT
C               CHARGE AND MAGNETIC FIELD ARE NEEDED
C               FOR CORRELATION TERMS (LAMBDA,Y),(LAMBDA,Z),(PHI,Y),(PHI,Z)
C               THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
C               TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
C               AND RS FOR FIXED X
C
C     IERR              NOT USED
C     SPX       SIGN OF X-COMPONENT OF PARTICLE MOMENTUM    INPUT
C
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL   PS,PC,H,RC,RS,CH,SPX
#endif
#include "geant321/trcom3.inc"
      DIMENSION PS(3),PC(3),H(3),RC(15),RS(15)
      DIMENSION UN(3),VN(3)
C
#if defined(CERNLIB_SINGLE)
      DATA CFACT8 / 2.997925 E-4 /
#endif
#if !defined(CERNLIB_SINGLE)
      DATA CFACT8 / 2.997925 D-4 /
#endif
C
      IERR=0
      PM=PS(1)
      TN(1)=1./SQRT(1.+PS(2)**2+PS(3)**2)
      IF(SPX.LT.0.) TN(1)=-TN(1)
      TN(2)=PS(2)*TN(1)
      TN(3)=PS(3)*TN(1)
C
      PC(1)=PS(1)
      PC(2)=ASIN(TN(3))
      IF (ABS (TN(1)) .LT. 1.E-30) TN(1) = 1.E-30
      PC(3) = ATAN2 (TN(2),TN(1))
C
      COSL=SQRT(ABS(1.-TN(3)**2))
      IF (COSL .LT. 1.E-30) COSL = 1.E-30
      COSL1=1./COSL
      UN(1)=-TN(2)*COSL1
      UN(2)=TN(1)*COSL1
C
      VN(1)=-TN(3)*UN(2)
      VN(2)=TN(3)*UN(1)
      VN(3)=COSL
C
      J=0
      DO 10 I=1,5
         DO 4 K=I,5
            J=J+1
            A(I,K)=0.
            A(K,I)=0.
            S(J)=RS(J)
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
      SINZ=-(H(1)*UN(1)+H(2)*UN(2)           )*HM
      COSZ= (H(1)*VN(1)+H(2)*VN(2)+H(3)*VN(3))*HM
      A(2,4)=-Q*TN(2)*SINZ
      A(2,5)=-Q*TN(3)*SINZ
      A(3,4)=-Q*TN(2)*COSZ*COSL1
      A(3,5)=-Q*TN(3)*COSZ*COSL1
C
    6 A(1,1)=1.
      A(2,2)=TN(1)*VN(2)
      A(2,3)=TN(1)*VN(3)
      A(3,2)=TN(1)*UN(2)*COSL1
      A(4,4)=UN(2)
      A(5,4)=VN(2)
      A(5,5)=VN(3)
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
*
* *** ERROR EXITS
*
      END
