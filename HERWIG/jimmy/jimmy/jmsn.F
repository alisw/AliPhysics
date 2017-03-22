      DOUBLE PRECISION FUNCTION JMSN( DUMMY1, DUMMY2 )
C --------------------------------------------------------------------
C  Purpose: To calculate the cross section for N and
C           only N scatters.
C --------------------------------------------------------------------

      IMPLICIT NONE

#include "jimmy.inc"


      INTEGER N, DUMMY2
      DOUBLE PRECISION JMSNFN, EPS, XSECN, DUMMY1
      INTEGER NDIM, MAXPTS, MINPTS, LENWRK, I
      PARAMETER(NDIM=2,MAXPTS=10000*NDIM,
     *   LENWRK=(NDIM+2)*(1+MAXPTS/(2**NDIM+2*NDIM*NDIM+2*NDIM+1)))

      REAL A(NDIM), B(NDIM), RESULT, ERROR

      EXTERNAL JMSNFN
      COMMON/BLOCK2/ XSECN, N       ! COMMON local to JMSNFN

      XSECN = DUMMY1
      N = DUMMY2
      DO 10 I=1,NDIM
        A(I)=0.0
        B(I)=1.0
  10  CONTINUE

      MINPTS = 0
      EPS    = 1.D-3
      FN_TYPE = 4
      CALL PARTN( NDIM, A, B, 1.0, 10000 )
      CALL INTGRL( NDIM, 0, 200, RESULT, ERROR )
      IF (JMBUG.GT.0) write(*,*) 'jmsn: result, error',result,error
      JMSN = DBLE(RESULT)

      RETURN
      END
