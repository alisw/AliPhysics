*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DCHPWS(N,C,A)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RCHPWS(N,C,A)
#endif
#include "gen/imp64.inc"
      CHARACTER NAMECP*(*),NAMEPC*(*)
      CHARACTER*80 ERRTXT
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAMECP = 'DCHPWS', NAMEPC = 'DPWCHS')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAMECP = 'RCHPWS', NAMEPC = 'RPWCHS')
#endif

      DIMENSION C(0:*),A(0:*),QU(0:101),QV(0:101)

      PARAMETER (R1 = 1, HF = R1/2)

      IF(N .LT. 0 .OR. N .GT. 100) THEN
       WRITE(ERRTXT,101) N
       CALL MTLPRT(NAMECP,'E408.1',ERRTXT)
      ELSEIF(N .EQ. 0) THEN
       A(0)=C(0)
      ELSE
       QU(0)=C(N)
       QU(1)=0
       QV(0)=C(N-1)
       QV(1)=2*QU(0)
       DO 2 K = 2,N
       TT=C(N-K)-QU(0)
       QU(K)=0
       DO 1 J = 0,K-1
       T=2*QV(J)-QU(J+1)
       QU(J)=QV(J)
       QV(J)=TT
       TT=T
    1  CONTINUE
       QV(K)=TT
    2  CONTINUE
      ENDIF
      A(0)=QV(0)
      DO 3 J = 1,N
      A(J)=QV(J)-QU(J-1)
    3 CONTINUE
      RETURN

#if defined(CERNLIB_DOUBLE)
      ENTRY DPWCHS(N,A,C)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY RPWCHS(N,A,C)
#endif

      IF(N .LT. 0 .OR. N .GT. 100) THEN
       WRITE(ERRTXT,101) N
       CALL MTLPRT(NAMEPC,'E408.1',ERRTXT)
      ELSEIF(N .EQ. 0) THEN
       QU(0)=A(0)
      ELSE
       QU(0)=A(N-1)
       QU(1)=A(N)
       IF(N .GE. 2) QU(2)=0
       DO 4 K = 2,N
       TT=QU(0)+HF*QU(2)
       QU(0)=A(N-K)+HF*QU(1)
       DO 5 J = 2,K
       T=HF*QU(J-1)
       IF(J .LE. K-2) T=T+HF*QU(J+1)
       QU(J-1)=TT
       TT=T
    5  CONTINUE
       QU(K)=T
    4  CONTINUE
      ENDIF
      DO 6 J = 0,N
      C(J)=QU(J)
    6 CONTINUE
      RETURN
  101 FORMAT('NUMBER OF TERMS  N = ',I5,'  < 0  OR  > 100')
      END
