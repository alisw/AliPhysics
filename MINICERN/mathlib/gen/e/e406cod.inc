*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
*
* e406cod.inc
*
      DIMENSION C(0:128),CC(0:128),AC(0:64),BC(0:64)

      PARAMETER (Z1 = 1, Z2 = 2, HF = Z1/2)

      ALFA=HF*(B-A)
      BETA=HF*(B+A)
      C1=F(A)
      C2=F(B)
      AC(0)=C2+C1
      AC(1)=C2-C1

      DO 1 I = 1,7
      I1=2**(I-1)
      I2=I1-1
      I3=2*I1
      C1=Z2/I1
      C2=PI/I1

      DO 2 J = 0,I2
    2 C(J)=F(ALFA*COS((J+HF)*C2)+BETA)

C    COMPUTE B-COEFFICIENTS

      DO 3 J = 0,I2
      F1=J*C2
      F2=-HF*F1

      C3=2*COS(F1)
      A2=0
      A1=0
      A0=C(I2)
      DO 4 K = I2-1,0,-1
      A2=A1
      A1=A0
    4 A0=C(K)+C3*A1-A2
    3 BC(J)=C1*(A0*COS(F1+F2)-A1*COS(F2))
      BC(I1)=0

C    COMPUTE NEW C-COEFFICIENTS

      DO 5 J = 0,I1
      C(J)=HF*(AC(J)+BC(J))
    5 C(I1+J)=HF*(AC(I1-J)-BC(I1-J))

C    CHECK IF CALCULATION IS FINISHED

      CMX=0
      DO 6 J = 0,I3
      CC(J)=ABS(C(J))
    6 CMX=MAX(CMX,CC(J))

      IF(CMX .GT. 0) THEN
       CMX=1/CMX
       CC(I3)=HF*CC(I3)
       A0=CC(I2)*CMX
       A1=CC(I1)*CMX
       DO 10 J = I1+2,I3
       A2=CC(J)*CMX
       IF(A0 .LE. EPS .AND. A1. LE. EPS .AND. A2 .LE. EPS) GO TO 9

       A0=A1
   10  A1=A2
      ENDIF

C    DOUBLE THE NUMBER OF COEFFICIENTS.

      IF(I .EQ. 7) GO TO 1

      DO 12 J = 0,I3
   12 AC(J)=C(J)

    1 CONTINUE

C    REQUIRED ACCURACY NOT OBTAINED

      NC=64
      DELTA=0
      DO 13 J = 60,NC
   13 DELTA=DELTA+ABS(C(J))
      C(0)=HF*C(0)
      CALL MTLPRT(NAME,'E406.1','REQUIRED ACCURACY NOT OBTAINED')
      RETURN

C    REQUIRED ACCURACY OBTAINED
C    SUM NEGLECTED TERMS IN EXPANSION

    9 DELTA=0
      DO 7 K = J,I3
    7 DELTA=DELTA+CC(K)

C    CHECK IF FURTHER REDUCTION OF COEFFICIENTS IS POSSIBLE.

      NC=J-1
      REST=EPS-DELTA
      IF(REST .GT. 0) THEN
    8  IF(CC(NC) .GE. REST) GO TO 14
       DELTA=DELTA+CC(NC)
       REST=REST-CC(NC)
       NC=NC-1
       GO TO 8
      ENDIF
   14 C(0)=HF*C(0)
      RETURN
      END
