*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DSMPLX(A,B,C,Z0,IDA,M,M1,N,N1,LW,IDW,W,X,Z,ITYPE)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RSMPLX(A,B,C,Z0,IDA,M,M1,N,N1,LW,IDW,W,X,Z,ITYPE)
#endif
#include "gen/imp64.inc"

      LOGICAL L1,L2,L3
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DSMPLX')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RSMPLX')
#endif
      PARAMETER (R10 = 10, IE0 = 15)

      DIMENSION A(IDA,*),B(*),C(*),W(*),X(*),LW(IDW,*)

      IF(M .LE. 0 .OR. N .LE. 0) THEN
       WRITE(ERRTXT,101) M,N
       CALL MTLPRT(NAME,'H101.1',ERRTXT)
       RETURN
      ENDIF
      IF(.NOT.(0 .LE. M1 .AND. M1 .LE. M .AND.
     1         0 .LE. N1 .AND. N1 .LE. N)) THEN
       WRITE(ERRTXT,102) M,N,M1,N1
       CALL MTLPRT(NAME,'H101.2',ERRTXT)
       RETURN
      ENDIF

      EPS=0
      DO 5 I = 1,M
      DO 5 K = 1,N
    5 EPS=EPS+ABS(A(I,K))
      EPSL=LOG10(2*EPS/(M*N))
      IEXP=INT(EPSL)-IE0
      IF(EPSL .LT. 0) IEXP=IEXP-1
      EPS=R10**IEXP

      Z=Z0
      DO 10 I = 1,M
      LW(I,1)=0
      IF(I .GT. M1) LW(I,1)=1
   10 LW(I,4)=I
      DO 11 K = 1,N
      LW(K,2)=0
      IF(K .GT. N1) LW(K,2)=-1
      LW(K,5)=M+K
   11 LW(K,3)=0

C     LW(I,1), LW(K,2) =  0   x >= 0
C     LW(I,1), LW(K,2) = -1   x = 0
C     LW(I,1), LW(K,2) =  1   x unrestricted
C     LW(I,1), LW(K,2) =  2   x linearly independent, no influence
C     LW(I,1)          = -2   x unrestricted, cannot be eliminated

C     Elimination of the equality constraints from the basis

      IF(N1 .NE. N) THEN
   55  KP=0
       DO 60 K = 1,N
       IF(LW(K,2) .EQ. -1) KP=K
   60  CONTINUE

C     If KP = 0, no more equality constraints in the basis

       IF(KP .NE. 0) THEN
        DMAX=0
        DO 61  I = 1, M
        IF(LW(I,1) .NE. -1) THEN
         DMAX=MAX(ABS(A(I,KP)),DMAX)
         IF(ABS(A(I,KP)) .EQ. DMAX) IP=I
        ENDIF
   61   CONTINUE
        IF(DMAX .LT. EPS) DMAX=0
        IF(DMAX .EQ. 0) THEN
         IF(ABS(C(KP)) .LT. EPS) C(KP)=0
         IF(C(KP) .EQ. 0) THEN
          LW(KP,2)=1
          GO TO 55
         ENDIF
         ITYPE=4
         GO TO 9
        ELSE

C     Homogeneous part of at least 2 equ. is linearly dependent
C     Inhomogeneous part is contradictory.

         CALL H101S1(A,B,C,Z,M,N,IDA,IP,KP,LW,IDW,EPS)
         LW(KP,2)=LW(IP,1)
         LW(IP,1)=-1
         GO TO 55
        ENDIF
       ENDIF

C     Non-basic variables only equation variables?

       IND=0
       DO 81 I = 1,M
       IF(LW(I,1) .NE. -1)  IND=IND+1
   81  CONTINUE

C     If non-basic variables only equation variables,
C     is constraint violated?

       IF(IND .EQ. 0) THEN
        DMIN=0
        DO 85  K = 1, N
        IF(LW(K,2) .NE. 1) THEN
         IF(ABS(C(K)) .LT. EPS) C(K)=0
         DMIN=MIN(C(K),DMIN)
        ENDIF
   85   CONTINUE
        ITYPE=1
        IF(DMIN .LT. 0) ITYPE=4
        GO TO 9
       ENDIF
      ENDIF

C     Eliminate unrestricted variables

      IF(M1 .EQ. M) GO TO 200
  111 IP=0
      DO 105 I = 1,M
      IF(LW(I,1) .EQ. 1) IP=I
  105 CONTINUE
      IF(IP .EQ. 0) GO TO 200

C     If there are no free variables, or if they have been
C     exchanged with equations, go to 200.
C     Present tableau represents feasable initial solution?

      DMAX=0
      DMIN=0
      DO 110  K = 1, N
      IF(LW(K,2) .NE. 1) THEN
       IF(ABS(C(K)) .LT. EPS) C(K)=0
       DMIN=MIN(C(K),DMIN)
       DMAX=MAX(C(K),DMAX)
      ENDIF
  110 CONTINUE
      IND=0
      IF(DMIN .GE. 0) IND=1
      Q=DMAX

C     If tableau represents an initial solution (IND = 1),
C     elimination of remaining unrestricted variables.

      DMAX=0
      DMIN=0
      DO 112  K = 1, N
      IF(LW(K,2) .NE. 1) THEN
       IF(ABS(A(IP,K)) .LT. EPS) A(IP,K)=0
       DMAX=MAX(A(IP,K),DMAX)
       IF(DMAX .EQ. A(IP,K)) KPMAX=K
       DMIN=MIN(A(IP,K),DMIN)
       IF(DMIN .EQ. A(IP,K)) KPMIN=K
      ENDIF
  112 CONTINUE
      IF(DMAX .NE. 0 .OR. DMIN .NE. 0) THEN
       IF(IND .NE. 0) THEN
        IF(ABS(B(IP)) .LT. EPS)  B(IP)=0
        IF(B(IP) .LT. 0 .AND. DMAX .GT. 0) THEN
         Q=Q/DMAX
         DO 150 K = 1,N
         IF(LW(K,2) .NE. 1 .AND. A(IP,K) .GT. 0) THEN
          QQ=C(K)/A(IP,K)
          Q=MIN(QQ,Q)
          IF(Q .EQ. QQ) KP=K
         ENDIF
  150    CONTINUE
        ELSEIF(B(IP) .GT. 0 .AND. DMIN .LT. 0) THEN
         Q=Q/DMIN
         DO 170 K = 1,N
         IF(LW(K,2) .NE. 1 .AND. A(IP,K) .LT. 0) THEN
          QQ=C(K)/A(IP,K)
          Q=MAX(QQ,Q)
          IF(Q .EQ. QQ) KP=K
         ENDIF
  170    CONTINUE
        ELSEIF(B(IP) .EQ. 0) THEN
         Q=Q/MAX(DMAX,ABS(DMIN))
         DO 190 K = 1,N
         IF(LW(K,2) .NE. 1 .AND. A(IP,K) .NE. 0) THEN
          QQ=ABS(C(K)/A(IP,K))
          Q=MIN(QQ,Q)
          IF(Q .EQ. QQ) KP=K
         ENDIF
  190    CONTINUE
        ENDIF
       ELSE
        KP=KPMAX
        IF(DMAX .LT. ABS(DMIN)) KP=KPMIN
       ENDIF
       CALL H101S1(A,B,C,Z,M,N,IDA,IP,KP,LW,IDW,EPS)
       LW(KP,2)=1
       LW(IP,1)=0
      ELSE

C     If all pivot elements for free variable = 0, it has to
C     remain out of the basis (LW(.,1) = -2).

       LW(IP,1)=-2
      ENDIF
      GO TO 111

C     Feasable initial solution by mutiphase method.

  200 IF(M1 .EQ. M .AND. N1 .EQ. N) GO TO 250
      KP=0
      DO 201 K = 1,N
      IF(LW(K,2) .EQ. 0) KP=KP+1
  201 CONTINUE
      IF(KP .EQ. 0) THEN
       IP=0
       DO 202 I = 1,M
       IF(LW(I,1) .NE. -1) IP=IP+1
  202  CONTINUE
       IF(IP .EQ. 0) THEN
        ITYPE=1
        GO TO 9
       ENDIF

       L1=.FALSE.
       L2=.FALSE.
       L3=.FALSE.
       DO 206 I = 1,M
       IF(LW(I,1) .NE. -1) THEN
        IF(ABS(B(I)) .LT. EPS) B(I)=0
        L1=LW(I,1) .EQ.  0 .AND. B(I) .GT. 0
        L2=LW(I,1) .EQ.  0 .AND. B(I) .EQ. 0
        L3=LW(I,1) .EQ.  0 .AND. B(I) .LT. 0
        L2=LW(I,1) .EQ. -2 .AND. B(I) .EQ. 0
        L3=LW(I,1) .EQ. -2 .AND. B(I) .NE. 0
       ENDIF
  206  CONTINUE
       IF(L1) ITYPE=1
       IF(L2) ITYPE=2
       IF(L3) ITYPE=3
       GO TO 9
      ENDIF

C     Only unrestricted variables in the basis,
C     hence final tableau obtained.

      IP=0
      DO 211 I = 1,M
      IF(LW(I,1) .EQ. 0)  IP=IP+1
  211 CONTINUE
      IF(IP .NE. 0) GO TO 250
      DMIN=0
      DO 212 K = 1,N
      IF(LW(K,2) .NE. 1) THEN
       IF(ABS(C(K)) .LT. EPS) C(K)=0
       DMIN=MIN(C(K),DMIN)
      ENDIF
  212 CONTINUE
      IF(DMIN .LT. 0) THEN
       ITYPE=4
       GO TO 9
      ENDIF

      IP=0
      DO 216 I = 1,M
      IF(LW(I,1) .EQ. -2) IP=IP+1
  216 CONTINUE
      IF(IP .EQ. 0) THEN
       ITYPE=1
       GO TO 9
      ENDIF

      L2=.FALSE.
      L3=.FALSE.
      DO 225 I = 1,M
      IF(LW(I,1) .NE. -1) THEN
       IF(ABS(B(I)) .LT. EPS) B(I)=0
       L2=B(I) .EQ. 0
       L3=.NOT.L2
      ENDIF
  225 CONTINUE
      IF(L2) ITYPE=2
      IF(L3) ITYPE=3
      GO TO 9

C     Variables out of the basis either unrestricted or 0,
C     they cannot be exchanged, hence final tableau obtained.
C
C     Tableau representing initial solution?

  250 DMIN=0
      DO 255 K = 1,N
      IF(LW(K,2) .NE. 1) THEN
       IF(ABS(C(K)) .LT. EPS) C(K)=0
       DMIN=MIN(C(K),DMIN)
       IF(C(K) .EQ. DMIN) KK=K
      ENDIF
  255 CONTINUE

C     All C(K) for constraints X >= 0, hence initial solution found.
C     Otherwise column with index KK new objective function.

      IF(DMIN .LT. 0) THEN
  270  DMIN=0
       DO 300 I = 1,M
       IF(LW(I,1) .EQ. 0) THEN
        IF(ABS(A(I,KK)) .LT. EPS) A(I,KK)=0
        DMIN=MIN(A(I,KK),DMIN)
        IF(A(I,KK) .EQ. DMIN) IP=I
       ENDIF
  300  CONTINUE
       IF(DMIN .GE. 0) THEN
        ITYPE=4
        GO TO 9
       ENDIF

       CALL H101S2(A,B,C,M,N,IDA,IP,KP,KK,LW,IDW,W,X,EPS,0)
       IF(KP .EQ. 0) KP=KK
       CALL H101S1(A,B,C,Z,M,N,IDA,IP,KP,LW,IDW,EPS)

       IF(ABS(C(KK)) .LT. EPS) C(KK)=0
       IF(C(KK) .LT. 0) GO TO 270
       GO TO 250
      ENDIF

C     Initial solution found, hence algorithm may stop.
C     Maximum obtained?

  500 KASE=0
      DO 510 I = 1,M
      IF(LW(I,1) .EQ. -2) THEN
       IF(ABS(B(I)) .LT. EPS) B(I)=0
       IF(B(I) .NE. 0) THEN
        ITYPE=3
        GO TO 9
       ENDIF
       KASE=1
      ENDIF
  510 CONTINUE

C     KASE = 1: There is a free variable out of the basis
C     which does not influence the maximum (ITYPE = 2).
C     Take out negative B(I).

      IP=0
      DO 555 I = 1,M
      IF(LW(I,1) .EQ. 0) THEN
       IF(ABS(B(I)) .LT. EPS) B(I)=0
       IF(B(I) .LT. 0)  IP=I
      ENDIF
  555 CONTINUE

C     All B(I) non-negative, hence final tableau obtained
C     Solution unique?

      IF(IP .NE. 0) THEN
       CALL H101S2(A,B,C,M,N,IDA,IP,KP,0,LW,IDW,W,X,EPS,0)
       IF(KP .EQ. 0) THEN
        ITYPE=3
        GO TO 9
       ENDIF
       CALL H101S1(A,B,C,Z,M,N,IDA,IP,KP,LW,IDW,EPS)
       GO TO 500
      ENDIF

      IF(KASE .NE. 0) THEN
       ITYPE=2
       GO TO 9
      ENDIF
      IND=0
      DO 562 K = 1,N
      IF(LW(K,2) .EQ. 0) THEN
       IF(ABS(C(K)) .LT. EPS) C(K)=0
       IF(C(K) .EQ. 0) THEN
        IND=IND+1
        LW(K,2)=2
       ENDIF
      ENDIF
  562 CONTINUE
      KASE=IND
      IND=0
      DO 565 I = 1,M
      IF(LW(I,1) .EQ. 0) THEN
      IF(ABS(B(I)) .LT. EPS)  B(I)=0
       IF(B(I) .EQ. 0) THEN
        IND=IND+1
        LW(I,1)=2
       ENDIF
      ENDIF
  565 CONTINUE
      IF(IND .EQ. 0) THEN
       ITYPE=1
       GO TO 9
      ELSEIF(KASE .EQ. 0) THEN
       ITYPE=2
       GO TO 9
      ENDIF

C     Some C(K) = 0 and some B(I) = 0 in the final tableau.
C     Solution unique?

  575 DO 577 I = 1,M
      IF(LW(I,1) .EQ. 2) THEN
       DMAX=0
       DO 576 K = 1, N
       IF(LW(K,2) .EQ. 2) THEN
        IF(ABS(A(I,K)) .LT. EPS) A(I,K)=0
        DMAX=MAX(A(I,K),DMAX)
       ENDIF
  576  CONTINUE
       IF(DMAX .LE. 0) THEN
        ITYPE=2
        GO TO 9
       ENDIF
      ENDIF
  577 CONTINUE
      DO 581 K = 1,N
      IF(LW(K,2) .EQ. 2) THEN
       DMIN=1
       DO 578 I = 1,M
       IF(LW(I,1) .EQ. 2) DMIN=MIN(A(I,K),DMIN)
  578  CONTINUE
       IF(DMIN .GT. 0) THEN
        ITYPE=1
        GO TO 9
       ELSEIF(DMIN .EQ. 0) THEN
        DO 580  I = 1, M
        IF(LW(I,1) .EQ. 2 .AND. A(I,K) .GT. 0) LW(I,1)=-2
  580   CONTINUE
        LW(K,2)=1
       ENDIF
      ENDIF
  581 CONTINUE
      NC=0
      DO 582 K = 1,N
      IF(LW(K,2) .EQ. 2) NC=NC+1
  582 CONTINUE
      NR=0
      DO 583 I = 1,M
      IF(LW(I,1) .EQ. 2) THEN
       NR=NR+1
       IP=I
      ENDIF
  583 CONTINUE
      IF(NR .EQ. 0 .OR. NC .EQ. 0) THEN
       ITYPE=1
       GO TO 9
      ENDIF

      CALL H101S2(A,B,C,M,N,IDA,IP,KP,0,LW,IDW,W,X,EPS,2)

      IF(KP .EQ. 0) THEN
       ITYPE=2
       GO TO 9
      ENDIF

      CALL H101S1(A,B,C,Z,M,N,IDA,IP,KP,LW,IDW,EPS)
      GO TO 575

    9 DO 1 I = 1,M
    1 X(LW(I,4))=0
      DO 2 K = 1,N
    2 X(LW(K,5))=C(K)
      RETURN
  101 FORMAT('M = ',I4,'  OR  N = ',I4,'  ILLEGAL')
  102 FORMAT('M = ',I4,',  N = ',I4,':  M1 = ',I4,'  OR  N1 = ',I4,
     1       '  ILLEGAL')
      END
