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
      SUBROUTINE TRSCSP(PC,RC,PS,RS,H,CH,IERR,SPX)
C
C *** TRANSFORMS ERROR MATRIX
C     FROM   SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
C      TO  SPLINE VARIABLES (1/P,Y',Z',Y,Z)
C
C     Authors: A. Haas and W. Wittek
C
C
C *** PC(3)     1/P,LAMBDA,PHI                          INPUT
C     PS(3)     1/P,Y',Z'                              OUTPUT
C     H(3)      MAGNETIC FIELD                          INPUT
C     RC(15)    ERROR MATRIX IN   SC   VARIABLES        INPUT     (TRIANGLE)
C     RS(15)    ERROR MATRIX IN SPLINE VARIABLES       OUTPUT     (TRIANGLE)
C     CH        CHARGE OF PARTICLE                      INPUT
C               CHARGE AND MAGNETIC FIELD ARE NEEDED
C               FOR CORRELATION TERMS (Y',YT),(Y',ZT),(Z',YT),(Z',ZT)
C               THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
C               TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
C               AND RS FOR FIXED X
C
C     IERR  =   1       PARTICLE MOVES PERPENDICULAR TO X-AXIS
C                      ( Y',Z' ARE NOT DEFINED )
C     SPX       SIGN OF X-COMPONENT OF PARTICLE MOMENTUM   OUTPUT
C
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL   PC,PS,H,RC,RS,CH,SPX
#endif
#include "geant321/trcom3.inc"
      DIMENSION PC(3),PS(3),H(3),RC(15),RS(15)
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
      PM=PC(1)
      COSL=COS(PC(2))
      SINP=SIN(PC(3))
      COSP=COS(PC(3))
C
      TN(1)=COSL*COSP
      SPX=1.
      IF(TN(1).LT.0.) SPX=-1.
      IF(TN(1).EQ.0.) GO TO 901
      TN(2)=COSL*SINP
      TN(3)=SIN(PC(2))
C
      T1R=1./TN(1)
      PS(1)=PC(1)
      PS(2)=SINP/COSP
      PS(3)=TN(3)*T1R
C
      UN(1)=-SINP
      UN(2)=COSP
C
      VN(1)=-TN(3)*UN(2)
      VN(2)=TN(3)*UN(1)
      VN(3)=COSL
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
      SINZ=-(H(1)*UN(1)+H(2)*UN(2)           )*HM
      COSZ= (H(1)*VN(1)+H(2)*VN(2)+H(3)*VN(3))*HM
 
      T3R=Q*T1R**3
      A(2,4)=-UN(1)*(VN(3)*COSZ           )*T3R
      A(2,5)=-VN(1)*(VN(3)*COSZ           )*T3R
      A(3,4)=UN(1)*(VN(2)*COSZ-UN(2)*SINZ)*T3R
      A(3,5)=VN(1)*(VN(2)*COSZ-UN(2)*SINZ)*T3R
C
    6 T2R=T1R**2
C
      A(1,1)=1.
      A(2,3)=VN(3)*COSL*T2R
      A(3,2)=UN(2)*T2R
      A(3,3)=-VN(2)*COSL*T2R
      A(4,4)=VN(3)*T1R
      A(5,4)=-VN(2)*T1R
      A(5,5)=UN(2)*T1R
C
      CALL SSMT5T(A,S,S)
C
      J=0
      DO 25 I=1,5
         DO 20 K=I,5
            J=J+1
            RS(J)=S(J)
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
