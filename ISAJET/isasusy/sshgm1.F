#include "isajet/pilot.h"
      SUBROUTINE SSHGM1(TAU,IFFF,RFFF)
C-----------------------------------------------------------------------
C
C     This subroutine uses the tau variable of the Higgs Hunters'
C     Guide.  Many other authors, including the paper cited in 
C     Higgs Hunters' Guide (PR. D. 38(11): 3481) and Collider Physics
C     by Barger and Phillips use the variable lambda
C          LAMBDA = ( MASS OF PARTICLE IN LOOP / MASS OF HIGGS )**2
C          TAU = 4.0 * LAMBDA
C     Note also that what is defined as the f function by different 
C     authors varies by a constant factor.  For example,
C          f(Barger and Phillips) = -2 * f(Higgs Hunters' Guide)
C
C     Bisset's FFF
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      DOUBLE PRECISION TAU,IFFF,RFFF,ETAP,ETAM,PI
C
      PI=3.1415926536
      IFFF=0.0
      RFFF=0.0
      IF(TAU.GE.1.0) THEN
        RFFF=(DASIN(DSQRT(1.0/TAU)))**2
        IFFF=0.0D0
      ELSE IF (TAU.LT.1.0) THEN
        ETAP=1.0D0+SQRT(1.0D0-TAU)
        ETAM=1.0D0-SQRT(1.0D0-TAU)
        RFFF=-((DLOG(ETAP/ETAM))**2-PI**2)/4.0D0
        IFFF=PI*DLOG(ETAP/ETAM)/2.0D0
      ENDIF
      RETURN 
      END
