*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RCSPLN(N,X,M,Y,NDIM,MODE,A,B,C,D)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DCSPLN(N,X,M,Y,NDIM,MODE,A,B,C,D)
#include "gen/imp64.inc"
#endif
      CHARACTER NAMEN*(*),NAMET*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEN = 'RCSPLN', NAMET = 'RCSPLT')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEN = 'DCSPLN', NAMET = 'DCSPLT')
#endif
      LOGICAL LNT

      DIMENSION X(0:*),Y(0:NDIM,*)
      DIMENSION A(NDIM,*),B(NDIM,*),C(NDIM,*),D(NDIM,*)

      PARAMETER (Z1 = 1, C1 = 3*Z1/2, C2 = Z1/3, C3 = 2*Z1/3)
      PARAMETER (C4 = Z1/2, C5 = Z1/6, C6 = 2*Z1/15, C7 = 7*Z1/60)

      LNT=.FALSE.
      GO TO 50

#if !defined(CERNLIB_DOUBLE)
      ENTRY RCSPNT(N,X,M,Y,NDIM,MODE,A,B,C,D)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DCSPNT(N,X,M,Y,NDIM,MODE,A,B,C,D)
#endif

      LNT=.TRUE.
   50 IF(N .LT. 2 .OR. M .LT. 1 .OR. NDIM .LT. N .OR.
     1   MODE .NE. 0. AND. MODE .NE. 1) THEN
       IF(N .LT. 2) THEN
        WRITE(ERRTXT,101) N
        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.1',ERRTXT)
        IF(     LNT) CALL MTLPRT(NAMET,'E211.1',ERRTXT)
       ENDIF
       IF(M .LT. 1) THEN
        WRITE(ERRTXT,102) M
        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.2',ERRTXT)
        IF(     LNT) CALL MTLPRT(NAMET,'E211.2',ERRTXT)
       ENDIF
       IF(NDIM .LT. N) THEN
        WRITE(ERRTXT,103) NDIM,N
        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.3',ERRTXT)
        IF(     LNT) CALL MTLPRT(NAMET,'E211.3',ERRTXT)
       ENDIF
       IF(MODE .NE. 0 .AND. MODE .NE. 1) THEN
        WRITE(ERRTXT,104) MODE
        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.4',ERRTXT)
        IF(     LNT) CALL MTLPRT(NAMET,'E211.4',ERRTXT)
       ENDIF
       RETURN
      ENDIF
      DO 1 I = 1,N
    1 D(I,1)=X(I)-X(I-1)
      DO 2 K = 1,M
      DO 2 I = 1,N
    2 A(I,K)=(Y(I,K)-Y(I-1,K))/D(I,1)
      IF(MODE .EQ. 1) THEN
       IF(N .EQ. 2) THEN
        T1=1/(D(1,1)+D(2,1))
        DO 3 K = 1,M
    3   C(2,K)=T1*(A(2,K)-A(1,K))
       ELSE
        DO 4 I = 2,N
    4   B(I,1)=1/(D(I,1)+D(I-1,1))
        DO 5 K = 1,M
    5   C(1,K)=0
        B(1,1)=1

        DO 6 I = 2,N-1
        T1=3*B(I,1)
        T2=B(I,1)*D(I-1,1)
        T3=1/(2+T2*B(I-1,1))
        B(I,1)=(T2-1)*T3
        DO 6 K = 1,M
    6   C(I,K)=(T1*(A(I,K)-A(I-1,K))-T2*C(I-1,K))*T3

        T1=3*B(N,1)
        T2=B(N,1)*D(N-1,1)
        T3=1/(3-T2*(1-B(N-1,1)))
        DO 7 K = 1,M
    7   C(N,K)=(T1*(A(N,K)-A(N-1,K))-T2*C(N-1,K))*T3
       END IF

       DO 8 I = N-1,2,-1
       T1=B(I,1)
       DO 8 K = 1,M
    8  C(I,K)=T1*C(I+1,K)+C(I,K)
       DO 9 K = 1,M
    9  C(1,K)=C(2,K)
       IF(.NOT.LNT) THEN
        DO 10 K = M,1,-1
        B(1,K)=A(1,K)-C(2,K)*D(1,1)
        D(1,K)=0
        DO 11 I = 2,N-1
        B(I,K)=A(I,K)-C2*(C(I+1,K)+2*C(I,K))*D(I,1)
   11   D(I,K)=(C(I+1,K)-C(I,K))/(3*D(I,1))
        B(N,K)=A(N,K)-C(N,K)*D(N,1)
   10   D(N,K)=0
        DO 12 K = 1,M
        DO 12 I = 1,N
   12   A(I,K)=Y(I-1,K)
       ENDIF
      ELSE
       IF(N .EQ. 2) THEN
        T1=C1/(D(1,1)+D(2,1))
        DO 23 K = 1,M
   23   C(2,K)=T1*(A(2,K)-A(1,K))
       ELSE
        DO 24 I = 2,N
   24   B(I,1)=1/(D(I,1)+D(I-1,1))
        DO 25 K = 1,M
   25   C(1,K)=0
        B(1,1)=0

        DO 26 I = 2,N-1
        T1=3*B(I,1)
        T2=B(I,1)*D(I-1,1)
        T3=1/(2+T2*B(I-1,1))
        B(I,1)=(T2-1)*T3
        DO 26 K = 1,M
   26   C(I,K)=(T1*(A(I,K)-A(I-1,K))-T2*C(I-1,K))*T3

        T1=3*B(N,1)
        T2=B(N,1)*D(N-1,1)
        T3=1/(2+T2*B(N-1,1))
        DO 27 K = 1,M
   27   C(N,K)=(T1*(A(N,K)-A(N-1,K))-T2*C(N-1,K))*T3
       END IF

       DO 28 I = N-1,2,-1
       T1=B(I,1)
       DO 28 K = 1,M
   28  C(I,K)=T1*C(I+1,K)+C(I,K)
       DO 29 K = 1,M
   29  C(1,K)=0
       IF(.NOT.LNT) THEN
        T1=C2*D(1,1)
        T2=C2/D(1,1)
        T3=C3*D(N,1)
        T4=C2/D(N,1)
        DO 30 K = M,1,-1
        B(1,K)=A(1,K)-T1*C(2,K)
        D(1,K)=T2*C(2,K)
        DO 31 I = 2,N-1
        B(I,K)=A(I,K)-C2*(C(I+1,K)+2*C(I,K))*D(I,1)
   31   D(I,K)=(C(I+1,K)-C(I,K))/(3*D(I,1))
        B(N,K)=A(N,K)-T3*C(N,K)
   30   D(N,K)=-T4*C(N,K)
        DO 32 K = 1,M
        DO 32 I = 1,N
   32   A(I,K)=Y(I-1,K)
       ENDIF
      ENDIF

      IF(LNT) THEN
       DO 41 K = 1,M
       T1=D(1,1)**2
       A(1,K)=C4*(Y(0,K)+Y(1,K)-C5*(C(1,K)+C(2,K))*T1)*D(1,1)
       B(1,K)=C2*(Y(0,K)+C4*Y(1,K)-(C6*C(1,K)+C7*C(2,K))*T1)*T1
       DO 42 I = 2,N-1
       T1=D(I,1)**2
       A(I,K)=A(I-1,K)+
     1        C4*(Y(I-1,K)+Y(I,K)-C5*(C(I,K)+C(I+1,K))*T1)*D(I,1)
   42  B(I,K)=B(I-1,K)+C2*(Y(I-1,K)+C4*Y(I,K)-
     1                (C6*C(I,K)+C7*C(I+1,K))*T1)*T1+A(I-1,K)*D(I,1)
       T1=D(N,1)**2
       A(N,K)=A(N-1,K)+
     1        C4*(Y(N-1,K)+Y(N,K)-C(N,K)*C5*(1+MODE)*T1)*D(N,1)
   41  B(N,K)=B(N-1,K)+C2*(Y(N-1,K)+C4*Y(N,K)-
     1        C(N,K)*(C6+MODE*C7)*T1)*T1+A(N-1,K)*D(N,1)
      ENDIF
      RETURN
  101 FORMAT('ILLEGAL N = ',I5,' < 2')
  102 FORMAT('ILLEGAL M = ',I5,' < 1')
  103 FORMAT('ILLEGAL NDIM = ',I5,' < ',I5,' = N')
  104 FORMAT('ILLEGAL MODE = ',I5)
      END
