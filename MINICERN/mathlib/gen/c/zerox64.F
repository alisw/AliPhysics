*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DZEROX(A0,B0,EPS,MAXF,F,MODE)
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RZEROX(A0,B0,EPS,MAXF,F,MODE)
#endif
C     Based on
C
C        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
C        Guaranteed Convergence for Finding a Zero of a Function,
C        ACM Trans. Math. Software 1 (1975) 330-345.
C
C        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
#include "gen/imp64.inc"
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DZEROX')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RZEROX')
#endif
      LOGICAL LMT

      DIMENSION IM1(2),IM2(2),LMT(2)

      PARAMETER (Z1 = 1, HALF = Z1/2)

      DATA IM1 /2,3/, IM2 /-1,3/

      IF(MODE .NE. 1 .AND. MODE .NE. 2) THEN
       C=0
       WRITE(ERRTXT,101) MODE
       CALL MTLPRT(NAME,'C200.1',ERRTXT)
       GO TO 99
      ENDIF
      FA=F(B0)
      FB=F(A0)
      IF(FA*FB .GT. 0) THEN
       C=0
       WRITE(ERRTXT,102) A0,B0
       CALL MTLPRT(NAME,'C200.2',ERRTXT)
       GO TO 99
      ENDIF
      ATL=ABS(EPS)
      B=A0
      A=B0
      LMT(2)=.TRUE.
      MF=2
    1 C=A
      FC=FA
    2 IE=0
    3 IF(ABS(FC) .LT. ABS(FB)) THEN
       IF(C .NE. A) THEN
        D=A
        FD=FA
       END IF
       A=B
       B=C
       C=A
       FA=FB
       FB=FC
       FC=FA
      END IF
      TOL=ATL*(1+ABS(C))
      H=HALF*(C+B)
      HB=H-B
      IF(ABS(HB) .GT. TOL) THEN
       IF(IE .GT. IM1(MODE)) THEN
        W=HB
       ELSE
        TOL=TOL*SIGN(Z1,HB)
        P=(B-A)*FB
        LMT(1)=IE .LE. 1
        IF(LMT(MODE)) THEN
         Q=FA-FB
         LMT(2)=.FALSE.
        ELSE
         FDB=(FD-FB)/(D-B)
         FDA=(FD-FA)/(D-A)
         P=FDA*P
         Q=FDB*FA-FDA*FB
        END IF
        IF(P .LT. 0) THEN
         P=-P
         Q=-Q
        END IF
        IF(IE .EQ. IM2(MODE)) P=P+P
        IF(P .EQ. 0 .OR. P .LE. Q*TOL) THEN
         W=TOL
        ELSEIF(P .LT. HB*Q) THEN
         W=P/Q
        ELSE
         W=HB
        END IF
       END IF
       D=A
       A=B
       FD=FA
       FA=FB
       B=B+W
       MF=MF+1
       IF(MF .GT. MAXF) THEN
        CALL MTLPRT(NAME,'C200.3','TOO MANY FUNCTION CALLS')
        GO TO 99
       ENDIF
       FB=F(B)
       IF(FB .EQ. 0 .OR. SIGN(Z1,FC) .EQ. SIGN(Z1,FB)) GO TO 1
       IF(W .EQ. HB) GO TO 2
       IE=IE+1
       GO TO 3
      END IF
#if defined(CERNLIB_DOUBLE)
      DZEROX=C
#endif
#if !defined(CERNLIB_DOUBLE)
      RZEROX=C
#endif
   99 CONTINUE
      RETURN
  101 FORMAT('MODE = ',I3,' ILLEGAL')
  102 FORMAT('F(A) AND F(B) HAVE THE SAME SIGN, A = ',1P,D15.8,
     1       ', B = ',D15.8)
      END
#if !defined(CERNLIB_DOUBLE)
      FUNCTION ZEROX(A0,B0,EPS,MAXF,F,MODE)
      ZEROX=RZEROX(A0,B0,EPS,MAXF,F,MODE)
      END
#endif
