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
      SUBROUTINE DSRTNT(K,N,A,B,C,U1,V1,RES,LRL)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RSRTNT(K,N,A,B,C,U1,V1,RES,LRL)
#endif
#include "gen/imp64.inc"
      LOGICAL LRL,LLL
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DSRTNT')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RSRTNT')
#endif
      PARAMETER (R1 = 1, HF = R1/2)

      DIMENSION BK(0:4,0:4),SGN(0:4)

      DATA BK(0,0) /1/, (BK(1,J),J=0,1) /1,1/, (BK(2,J),J=0,2) /1,2,1/
      DATA (BK(3,J),J=0,3) /1,3,3,1/, (BK(4,J),J=0,4) /1,4,6,4,1/
      DATA (SGN(J),J=0,4) /1,-1,1,-1,1/

      P(X)=A+B*X+C*X**2
      RT(X)=SQRT(P(X))

      H=0
      LLL=ABS(A)+ABS(B)+ABS(C) .GT. 0
      IF(.NOT.LLL) THEN
       CALL MTLPRT(NAME,'B300.1','A = B = C = 0')
       GO TO 9
      ENDIF
      LLL=ABS(K) .LE. 3 .AND. (N .EQ. 1 .OR. N .EQ. 3)
      IF(.NOT.LLL) THEN
       WRITE(ERRTXT,102) K,N
       CALL MTLPRT(NAME,'B300.2',ERRTXT)
       GO TO 9
      ENDIF

      U=MIN(U1,V1)
      V=MAX(U1,V1)
      LLL=U .EQ. V
      IF(LLL) GO TO 9

      LLL=K .GE. 0 .OR. K .EQ. -1 .AND. U*V .NE. 0
     1             .OR. K .LE. -2 .AND. U*V .GT. 0
      IF(.NOT.LLL) GO TO 9

      DELTA=4*A*C-B**2
      RTD=SQRT(ABS(DELTA))
      IF(C .EQ. 0) THEN
       IF(B .LT. 0) THEN
        X0=-A/B
        LLL=U .LE. X0 .AND. V .LE. X0
       ELSEIF(B .EQ. 0) THEN
        LLL=A .GT. 0
       ELSE
        X0=-A/B
        LLL=U .GE. X0 .AND. V .GE. X0
       ENDIF
      ELSE
       IF(DELTA .GT. 0) THEN
        LLL=C .GT. 0
       ELSEIF(DELTA .EQ. 0) THEN
        IF(C .LT. 0) THEN
         LLL=.FALSE.
        ELSE
         X0=-B/(2*C)
         LLL=U .LT. X0 .AND. V .LT. X0 .OR. U .GT. X0 .AND. V .GT. X0
        ENDIF
       ELSE
        A1=(-B+RTD)/(2*C)
        A2=(-B-RTD)/(2*C)
        W1=MIN(A1,A2)
        W2=MAX(A1,A2)
        IF(C .LT. 0) THEN
         LLL=W1 .LE. U .AND. U .LE. W2 .AND. W1 .LE. V .AND. V .LE. W2
        ELSE
         LLL=U .LE. W1 .AND. V .LE. W1 .OR. U .GE. W2 .AND. V .GE. W2
        ENDIF
       ENDIF
      ENDIF
      IF(.NOT.LLL) GO TO 9

      IF(K .GE. 0) THEN
       IF(C .EQ. 0) THEN
        IF(B .EQ. 0) THEN
         K1=K+1
         H=(V**K1-U**K1)/(K1*SQRT(A)**N)
        ELSE
         N1=2-N
         XV=A+B*V
         XU=A+B*U
         IF(A .EQ. 0) THEN
          K1=2*K+N1
          H=2*(SQRT(XV)**K1-SQRT(XU)**K1)/(K1*B**(K+1))
         ELSE
          HV=SQRT(XV)**N1
          HU=SQRT(XU)**N1
          H1=-XV/A
          H2=-XU/A
          S=N1*(HV-HU)
          DO 1 J = 1,K
    1     S=S+BK(K,J)*(H1**J*HV-H2**J*HU)/(2*J+N1)
          H=2*(-A/B)**K*S/B
         ENDIF
        ENDIF
       ELSE
        IF(N .EQ. 1) THEN
         ASSIGN 11 TO JMP1
         GO TO 10
   11    IF(K .EQ. 0) THEN
          H=H
         ELSEIF(K .EQ. 1) THEN
          H=(RT(V)-RT(U)-HF*B*H)/C
         ELSEIF(K .EQ. 2) THEN
          H1=4*C
          H2=6*B
          H=((H1*V-H2)*RT(V)-(H1*U-H2)*RT(U)+(2*B**2-DELTA)*H)/(8*C**2)
         ELSEIF(K .EQ. 3) THEN
          H1=8*C**2
          H2=10*B*C
          G1=15*B**2
          G2=A*C
          H3=G1-16*G2
          H=(((H1*V-H2)*V+H3)*RT(V)-((H1*U-H2)*U+H3)*RT(U)-
     1       (HF*G1-18*G2)*B*H)/(24*C**3)
         ENDIF
        ELSE
         IF(DELTA .EQ. 0) THEN
          IF(B .EQ. 0) THEN
           IF(K .EQ. 2) THEN
            H=LOG(ABS(V/U))/SQRT(C)**3
           ELSE
            K1=K-2
            H=(V**K1-U**K1)/(K1*SQRT(C)**3)
           ENDIF
           H=SIGN(R1,U)*H
          ELSE
           X0=B/(2*C)
           XV=V+X0
           XU=U+X0
           HV=1/XV**2
           HU=1/XU**2
           H1=-XV/X0
           H2=-XU/X0
           S=HF*(HU-HV)
           DO 2 J = 1,K
           IF(J .NE. 2) THEN
            S=S+BK(K,J)*(H1**J*HV-H2**J*HU)/(J-2)
           ELSE
            S=S+BK(K,2)*LOG(ABS(XV/XU))/X0**2
           ENDIF
    2      CONTINUE
           H=(-X0)**K*S/SQRT(C)**3
           IF(U .LT. -X0) H=-H
          ENDIF
         ELSE
          IF(K .EQ. 0) THEN
           H1=2*C
           H=2*((H1*V+B)/RT(V)-(H1*U+B)/RT(U))/DELTA
          ELSEIF(K .EQ. 1) THEN
           H1=2*A
           H=2*((H1+B*U)/RT(U)-(H1+B*V)/RT(V))/DELTA
           LB1=11
          ELSEIF(K .EQ. 2) THEN
           ASSIGN 12 TO JMP1
           GO TO 10
   12      H1=DELTA-B**2
           H2=2*A*B
           H=(((H1*U-H2)/RT(U)-(H1*V-H2)/RT(V))/DELTA+H)/C
          ELSEIF(K .EQ. 3) THEN
           ASSIGN 13 TO JMP1
           GO TO 10
   13      H1=C*DELTA
           G1=A*C
           G2=3*B**2
           H2=B*(10*G1-G2)
           H3=A*(8*G1-G2)
           H=(2*(((H1*V+H2)*V+H3)/RT(V)-((H1*U+H2)*U+H3)/RT(U))/
     1       DELTA-3*B*H)/(2*C**2)
          ENDIF
         ENDIF
        ENDIF
       ENDIF
      ELSE
       IF(A .EQ. 0) THEN
        IF(B .EQ. 0) THEN
         K1=K-N+1
         H=SIGN(R1,U)*(V**K1-U**K1)/(K1*SQRT(C)**N)
        ELSE
         IF(C .EQ. 0) THEN
          K1=2*K-N+2
          H=2*(SQRT(B*V)**K1-SQRT(B*U)**K1)/(K1*B**(K+1))
         ELSE
          XV=SQRT(C+B/V)
          XU=SQRT(C+B/U)
          N1=N-2
          K1=-K+N1
          S=0
          DO 4 J = 0,K1
          KJ=2*J-N1
    4     S=S+SGN(J)*BK(K1,J)*(XU**KJ-XV**KJ)/(KJ*C**J)
          H=2*(-C/B)**K1*S/B
          IF(U .LT. 0 .AND. V .LT. 0) H=-H
         ENDIF
        ENDIF
       ELSE
        IF(N .EQ. 1) THEN
         ASSIGN 21 TO JMP2
         GO TO 20
   21    IF(K .EQ. -1) THEN
          H=H
         ELSEIF(K .EQ. -2) THEN
          H=(RT(U)/U-RT(V)/V-HF*B*H)/A
         ELSEIF(K .EQ. -3) THEN
          H1=6*B
          H2=4*A
          H=((H1*V-H2)*RT(V)/V**2-(H1*U-H2)*RT(U)/U**2+
     1       (3*B**2-H2*C)*H)/(8*A**2)
         ENDIF
        ELSE
         IF(DELTA .EQ. 0) THEN
          IF(C .EQ. 0) THEN
           IF(K .EQ. -1) THEN
            H=LOG(ABS(V/U))/SQRT(A)**3
           ELSE
            K1=K+1
            H=(V**K1-U**K1)/(K1*SQRT(A)**3)
           ENDIF
          ELSE
           X0=B/(2*C)
           XV=1+X0/V
           XU=1+X0/U
           K1=-K+1
           K2=-K-1
           S=0
           DO 3 J = 0,K1
           KJ=K2-J
           IF(KJ .NE. 0) THEN
            S=S+SGN(J)*BK(K1,J)*(XV**KJ-XU**KJ)/KJ
           ELSE
            S=S+SGN(K2)*BK(K1,K2)*LOG(ABS(XV/XU))
           ENDIF
    3      CONTINUE
           H=-S/(SQRT(C)**3*X0**(K1+1))
           IF(U. LT. -X0) H=-H
          ENDIF
         ELSE
          ASSIGN 22 TO JMP2
          GO TO 20
   22     IF(K .EQ. -1) THEN
           H1=B*C
           H2=B**2-2*A*C
           H=(2*((H1*U+H2)/RT(U)-(H1*V+H2)/RT(V))/DELTA+H)/A
          ELSEIF(K .EQ. -2) THEN
           G1=3*B**2
           G2=A*C
           H1=(G1-8*G2)*C
           H2=(G1-10*G2)*B
           H3=A*DELTA
           H=(((H1*V+H2-H3/V)/RT(V)-(H1*U+H2-H3/U)/RT(U))/DELTA
     1       -3*HF*B*H)/A**2
          ELSEIF(K .EQ. -3) THEN
           G1=A*DELTA
           G2=A*C
           G3=B**2
           G4=15*G3
           H1=2*A*G1
           H2=5*B*G1
           H3=(G4-62*G2)*G3+24*G2**2
           H4=B*C*(G4-52*G2)
           H=((((H2-H1/V)/V-H3-H4*V)/RT(V)-((H2-H1/U)/U-H3-H4*U)/RT(U))
     1        /DELTA+HF*(G4-12*G2)*H)/(4*A**3)
          ENDIF
         ENDIF
        ENDIF
       ENDIF
      ENDIF
      GO TO 9

   10 C2=2*C
      IF(DELTA .GT. 0 .OR. DELTA .LT. 0 .AND. C .GT. 0) THEN
       H=LOG(ABS((2*SQRT(C*P(V))+C2*V+B)/
     1           (2*SQRT(C*P(U))+C2*U+B)))/SQRT(C)
      ELSEIF(DELTA .EQ. 0) THEN
       H=ABS(LOG(ABS((C2*V+B)/(C2*U+B))))/SQRT(C)
      ELSE
       H=(ASIN((C2*U+B)/RTD)-ASIN((C2*V+B)/RTD))/SQRT(-C)
      ENDIF
      GO TO JMP1, (11,12,13)

   20 IF(C .EQ. 0) THEN
       IF(B .EQ. 0) THEN
        H=LOG(ABS(V/U))/SQRT(A)
       ELSE
        IF(A .LT. 0) THEN
         H=2*(ATAN(SQRT(-(A+B*V)/A))-ATAN(SQRT(-(A+B*U)/A)))/SQRT(-A)
        ELSE
         WA=SQRT(A)
         WU=SQRT(A+B*U)
         WV=SQRT(A+B*V)
         H=LOG(ABS((WV-WA)*(WU+WA)/((WV+WA)*(WU-WA))))/WA
        ENDIF
       ENDIF
      ELSE
       A2=2*A
       IF(DELTA .GT. 0 .OR. DELTA .LT. 0 .AND. A .GT. 0) THEN
        H=LOG(ABS((-2*SQRT(A*P(V))+B*V+A2)*U/
     1           ((-2*SQRT(A*P(U))+B*U+A2)*V)))/SQRT(A)
       ELSEIF(DELTA .EQ. 0) THEN
        H=LOG(ABS((B*U+A2)*V/((B*V+A2)*U)))/SQRT(A)
        IF(U*V .GT. 0) H=SIGN(H,U)
       ELSE
        H=(ASIN((B*V+A2)/(V*RTD))-ASIN((B*U+A2)/(U*RTD)))/SQRT(-A)
        IF(U .LT. 0 .AND. V .LT. 0) H=-H
       ENDIF
      ENDIF
      GO TO JMP2, (21,22)

    9 RES=SIGN(R1,V1-U1)*H
      LRL=LLL
      RETURN
  102 FORMAT('ILLEGAL VALUE(S) K = ',I5,', N = ',I5)
      END
