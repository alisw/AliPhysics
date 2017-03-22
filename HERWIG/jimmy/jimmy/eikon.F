      DOUBLE PRECISION FUNCTION EIKON( DUMMY )
C     --------------------------------------------------------------------
C     Purpose: To eikonalise the total hadronic cross section pass the
C     cross section to be eikonalised (XSECN) and the routine returns
C     the eikonalised version (EIKON)
C     --------------------------------------------------------------------
      IMPLICIT NONE

#include "jimmy.inc"

      DOUBLE PRECISION EPS, XSECN, DUMMY
      INTEGER NDIM, MAXPTS, MINPTS, LENWRK, I
      PARAMETER(NDIM=2,MAXPTS=10000*NDIM,
     & LENWRK=(NDIM+2)*(1+MAXPTS/(2**NDIM+2*NDIM*NDIM+2*NDIM+1)))
      REAL A(NDIM), B(NDIM), RESULT, ERROR

      COMMON/BLOCK1/XSECN               !COMMON local to eikfun

      XSECN = DUMMY
      DO 10 I=1,NDIM
        A(I)=0.0
        B(I)=1.0
  10  CONTINUE

      MINPTS = 0
      EPS    = 1.D-3
      FN_TYPE = 3
      CALL PARTN( NDIM, A, B, 1.0, 10000 )
      CALL INTGRL( NDIM, 0, 200, RESULT, ERROR )
      EIKON = DBLE(RESULT)

      RETURN
      END


