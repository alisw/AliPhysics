*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DBINOM(X,K)
#include "gen/imp64.inc"
C     PARAMETER (IMX = (2**31)-1)
      PARAMETER (IMX = 2147483647)
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION BINOM(X,K)
C     PARAMETER (IMX = (2**31)-1)
      PARAMETER (IMX = 70368744177663)
#endif
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT

      PARAMETER (NAME = 'KBINOM')

      IF(K .LT. 0) THEN
       H=0
      ELSEIF(K .EQ. 0) THEN
       H=1
      ELSE
       H=X
       A=H
       DO 1 I = 2,K
       A=A-1
    1  H=H*(A/I)
      ENDIF
#if defined(CERNLIB_DOUBLE)
      DBINOM=H
#endif
#if !defined(CERNLIB_DOUBLE)
      BINOM=H
#endif
      RETURN

      ENTRY KBINOM(N,K)

    3 IF(N .GE. 0) THEN
       N1=N
       IS=1
      ELSE
       N1=K-N-1
       IS=(-1)**K
      ENDIF
      IF(K .LT. 0 .OR. K .GT. N1) THEN
       IH=0
      ELSEIF(K .EQ. 0 .OR. K .EQ. N1) THEN
       IH=IS
      ELSE
       IF(K+K .LE. N1) THEN
        K1=K
       ELSE
        K1=N1-K
       ENDIF
       H=N1
       A=H
       DO 2 I = 2,K1
       A=A-1
    2  H=H*(A/I)
       IF(H .LE. IMX) THEN
        IH=IS*NINT(H)
       ELSE
        IH=0
        WRITE(ERRTXT,101) N,K
        CALL MTLPRT(NAME,'B100.1',ERRTXT)
       ENDIF
      ENDIF
      KBINOM=IH
      RETURN
  101 FORMAT('INTEGER RANGE EXCEEDED, N =',I10,',  K = ',I10)
      END

