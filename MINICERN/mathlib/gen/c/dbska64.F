*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:08  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DBSKA(X,IA,JA,NL,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE BSKA(X,IA,JA,NL,B)
#include "gen/imp64.inc"
#endif
      LOGICAL LEX
 
      CHARACTER NAME*(*),ENAM*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'BSKA/DBSKA', ENAM = 'EBSKA/DEBKA')
 
      PARAMETER (Z1 = 1, Z2 = 2, Z3 = 3, Z4 = 4)
      PARAMETER (Z12 = Z1/Z2, Z13 = Z1/Z3, Z14 = Z1/Z4, Z23 = Z2/Z3)
      PARAMETER (Z34 = Z3/Z4)
 
      DIMENSION B(0:*)
 
      PARAMETER (PI = 3.14159 26535 89793D0, PIV = PI/4)
 
      LEX=.FALSE.
      GO TO 9
 
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBKA(X,IA,JA,NL,B)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY EBKA(X,IA,JA,NL,B)
#endif
      LEX=.TRUE.
 
    9 MODE=10*IA+JA
      N=NL-1
      U=2/X
      IF(LEX) THEN
       IF(X .LE. 0) THEN
        N=0
        WRITE(ERRTXT,101) X
        CALL MTLPRT(ENAM,'C341.1',ERRTXT)
       ELSEIF(NL .LT. 0 .OR. NL .GT. 100) THEN
        N=0
        WRITE(ERRTXT,103) NL
        CALL MTLPRT(ENAM,'C341.3',ERRTXT)
       ELSEIF(IA .EQ. 0) THEN
        A=0
        B(0)=DEBSK0(X)
        B(1)=DEBSK1(X)
       ELSEIF(MODE .EQ. 12) THEN
        A=Z12
        B(0)=SQRT(PIV*U)
        B(1)=B(0)*(1+A*U)
       ELSEIF(MODE .EQ. 13) THEN
        A=Z13
        B(0)=DEBKR3(X,1)
        B(1)=DEBKR3(X,2)+A*U*B(0)
       ELSEIF(MODE .EQ. 14) THEN
        A=Z14
        B(0)=DEBKR4(X,1)
        B(1)=DEBKR4(X,3)+A*U*B(0)
       ELSEIF(MODE .EQ. 23) THEN
        A=Z23
        B(0)=DEBKR3(X,2)
        B(1)=DEBKR3(X,1)+A*U*B(0)
       ELSEIF(MODE .EQ. 34) THEN
        A=Z34
        B(0)=DEBKR4(X,3)
        B(1)=DEBKR4(X,1)+A*U*B(0)
       ELSE
        N=0
        WRITE(ERRTXT,102) IA,JA
        CALL MTLPRT(ENAM,'C341.2',ERRTXT)
       ENDIF
      ELSE
       IF(X .LE. 0) THEN
        N=0
        WRITE(ERRTXT,101) X
        CALL MTLPRT(NAME,'C341.1',ERRTXT)
       ELSEIF(NL .LT. 0 .OR. NL .GT. 100) THEN
        N=0
        WRITE(ERRTXT,103) NL
        CALL MTLPRT(NAME,'C341.3',ERRTXT)
       ELSEIF(IA .EQ. 0) THEN
        A=0
        B(0)=DBESK0(X)
        B(1)=DBESK1(X)
       ELSEIF(MODE .EQ. 12) THEN
        A=Z12
        B(0)=EXP(-X)*SQRT(PIV*U)
        B(1)=B(0)*(1+A*U)
       ELSEIF(MODE .EQ. 13) THEN
        A=Z13
        B(0)=DBSKR3(X,1)
        B(1)=DBSKR3(X,2)+A*U*B(0)
       ELSEIF(MODE .EQ. 14) THEN
        A=Z14
        B(0)=DBSKR4(X,1)
        B(1)=DBSKR4(X,3)+A*U*B(0)
       ELSEIF(MODE .EQ. 23) THEN
        A=Z23
        B(0)=DBSKR3(X,2)
        B(1)=DBSKR3(X,1)+A*U*B(0)
       ELSEIF(MODE .EQ. 34) THEN
        A=Z34
        B(0)=DBSKR4(X,3)
        B(1)=DBSKR4(X,1)+A*U*B(0)
       ELSE
        N=0
        WRITE(ERRTXT,102) IA,JA
        CALL MTLPRT(NAME,'C341.2',ERRTXT)
       ENDIF
      ENDIF
      DO 1 J = 1,N
      A=A+1
    1 B(J+1)=B(J-1)+A*U*B(J)
      RETURN
  101 FORMAT('NON-POSITIVE ARGUMENT  X = ',E15.6)
  102 FORMAT('PAIR (IA,JA) = (',I5,I5,')  ILLEGAL')
  103 FORMAT('ILLEGAL  NL =',I5)
      END
