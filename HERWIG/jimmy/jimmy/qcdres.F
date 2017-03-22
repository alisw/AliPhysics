      DOUBLE PRECISION FUNCTION QCDRES( Z )
C     ---------------------------------------------------------------------
C     Purpose: To calculate the inclusive hadron(or photon)+hadron(or
C     photon) => JETS 
C     ---------------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"

      INTEGER I, MINPTS, NDIM, MAXPTS, LENWRK

C     Set parameters for integration
      PARAMETER(NDIM=3,MAXPTS=10000*NDIM,
     *   LENWRK=(NDIM+2)*(1+MAXPTS/(2**NDIM+2*NDIM*NDIM+2*NDIM+1)))
      DOUBLE PRECISION Z, EPS
      EXTERNAL JMKERN
      REAL A(NDIM), B(NDIM)
      REAL RESULT, ERROR

C     INITIALISE INTEGRATION LIMITS (ALL 0 TO 1)
      DO 10 I=1, NDIM
        A(I)=0.0
        B(I)=1.0
  10  CONTINUE

      MINPTS=0
      EPS=1.0D-3
C     Store the current z value
      YGAMMA=Z

      FN_TYPE=101
      CALL PARTN( NDIM, A, B, 1.0, 10000 )
      CALL INTGRL( NDIM, 0, 200, RESULT, ERROR )
      QCDRES=DBLE(RESULT)

C     Re-include the PHAD factor (JMKERN divides out the probability of the
C     photon becoming a hadronic = 1/phad)
      QCDRES = QCDRES/PHAD

      RETURN

      END







