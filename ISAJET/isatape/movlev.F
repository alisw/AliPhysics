#include "isajet/pilot.h"
#if defined(CERNLIB_MOVEFTN)
      SUBROUTINE MOVLEV(A,B,N)
C
C          Replacement for CDC system routine.
C          Move N consecutive locations from A to B.
C          Ver. 7.02: Separate entry points for real, integer, logical
C                     to comply with strict Fortran standard.
C                     Hence incompatible with CDC -- CDC version is 
C                     therefore obsolete, like the computer.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
      INTEGER N,I
      REAL A(N),B(N)
      INTEGER IA(N),IB(N)
      LOGICAL LA(N),LB(N)
C          Reals
      DO 100 I=1,N
        B(I)=A(I)
100   CONTINUE
      RETURN
C          Integers
      ENTRY MOVLEI(IA,IB,N)
      DO 200 I=1,N
        IB(I)=IA(I)
200   CONTINUE
      RETURN
C          Logicals
      ENTRY MOVLEL(LA,LB,N)
      DO 300 I=1,N
        LB(I)=LA(I)
300   CONTINUE
      RETURN
      END
#endif
