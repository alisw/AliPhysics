*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE WPOLYZ(C,N,MAXIT,Z0,H)
C
#include "gen/imp64.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'WPOLYZ')
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE CPOLYZ(C,N,MAXIT,Z0,H)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'CPOLYZ')
#endif
C
C     Based on
C      T. Pomentale, Homotopy methods for polynomial equations,
C      J. Inst. Maths Applics 13 (1974) 201-213

C     Revised version of a routine originally written by T. Pomentale

#include "gen/defc64.inc"
     + C(0:*),Z0(*),F(100),B(100),Z,CN,CD,FZ,FF,F1,F2,V,ZZ,Z1,ZV,IU
      LOGICAL LMD
      DIMENSION H(*)
      CHARACTER*80 ERRTXT

      PARAMETER (EPS = 1D-14)
#if defined(CERNLIB_CRAY)
      PARAMETER (IU = (0E0,1E0))
#endif
#if !defined(CERNLIB_CRAY)
      PARAMETER (IU = (0D0,1D0))
#endif
      PARAMETER (R0 = 0, R1 = 1, BIG1 = 1D7, BIG2 = BIG1/10, HF = R1/2)
      PARAMETER (FAC1 = 9*R1/10, FAC2 = 2+100*FAC1*EPS)
      PARAMETER (PI = 3.14159 26535 89793 24D0, PI2 = 2*PI)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
#if defined(CERNLIB_DOUBLE)
      TST(Z1,ZV)=DREAL(Z1)*DREAL(ZV)+DIMAG(Z1)*DIMAG(ZV)
#endif
#if !defined(CERNLIB_DOUBLE)
      TST(Z1,ZV)=REAL(Z1)*REAL(ZV)+AIMAG(Z1)*AIMAG(ZV)
#endif

      IF(N .LE. 0) RETURN
      IF(C(0) .EQ. 0) THEN
       CALL MTLPRT(NAME,'C209.1','A(0) = 0')
       RETURN
      ENDIF
      IF(N .EQ. 1) THEN
       Z0(1)=-C(1)/C(0)
       H(1)=0
       RETURN
      ENDIF
      IF(N .EQ. 2) THEN
       H(1)=0
       H(2)=0
       Z0(2)=(-C(1)+SQRT(C(1)**2-4*C(0)*C(2)))/(2*C(0))
       Z0(1)=-Z0(2)-C(1)/C(0)
       RETURN
      ENDIF
      LMD=.FALSE.
      DO 1 K = 1,N
    1 LMD=LMD .OR. Z0(K) .NE. 0
      M=0
      WN=SQRT(N+R0)
      Z=1+IU
      N1=N

      DO 3 K = 1,N+1
      H(K)=1
    3 B(K)=C(K-1)
    6 DO 10 K = 1,N-2
      IF(B(N1+1) .EQ. 0) GO TO 18
      IF(LMD) Z=Z0(K)
      IF(M .GE. MAXIT) THEN
       CALL MTLPRT(NAME,'C209.2','TOO MANY ITERATIONS')
       RETURN
      ENDIF

      S=BIG1
      DO 61 I = 1,N1+1
   61 F(I)=B(I)
      ZA=ABS(Z)
      HH=HF*ABS(B(1))
      DO 62 I = 2,N1+1
      F(I)=Z*F(I-1)+F(I)
      D=ABS(F(I))
   62 HH=ZA*HH+D
      DO 63 L = 1,2
      DO 63 I = 2,N1-L+1
   63 F(I)=Z*F(I-1)+F(I)
      M=M+1
      IF(D .LE. EPS*(FAC2*HH-D)) GO TO 13

    7 D1=S
    8 NI=0
      V=Z
      FZ=F(N1+1)
      FF=FZ
      F1=F(N1)
      F2=F(N1-1)

    9 NI=NI+1
      IH=0
      Z=V
      F(N1+1)=FF
      F(N1)=F1
      F(N1-1)=F2
   11 IH=IH+1
      IF(IH .GT. NI) GO TO 9
      CN=F(N1+1)-(NI-IH)*FZ/(NI+R0)
      CD=F(N1)**2-2*F(N1-1)*CN
      CN=CN*F(N1)
      IF(ABS(CN) .GT. 10*ABS(CD)) GO TO 9
      Z=Z-CN/CD
      IF(M .GE. MAXIT) THEN
       CALL MTLPRT(NAME,'C209.2','TOO MANY ITERATIONS')
       RETURN
      ENDIF

      DO 71 I = 1,N1+1
   71 F(I)=B(I)
      ZA=ABS(Z)
      HH=HF*ABS(B(1))
      DO 72 I = 2,N1+1
      F(I)=Z*F(I-1)+F(I)
      D=ABS(F(I))
   72 HH=ZA*HH+D
      DO 73 L = 1,2
      DO 73 I = 2,N1-L+1
   73 F(I)=Z*F(I-1)+F(I)
      M=M+1
      ERR=EPS*(FAC2*HH-D)
      IF(D .LE. ERR) GO TO 13
      S=ABS(F(N1+1))
      IF(S .LE. FAC1*D1) GO TO 7
      IF(S .LT. BIG2*ERR) GO TO 13
      S1=ABS(F(N1))
      S2=ABS(F(N1-1))
      S3=1+ABS(Z)
      IF(2*S1 .GT. S2*S3 .OR. 10*S .LT. S1*S3) GO TO 11

      ZZ=Z
      KK=2
      M1=1
      S=S1/S2
   40 DT=PI2/KK
      DO 42 J=1,KK,M1
      ZV=S*EXP(IU*(J*DT))
      Z=ZZ+ZV
      IF(M .GE. MAXIT) GO TO 44
      DO 41 I = 1,N1+1
   41 F(I)=B(I)
      DO 43 L = 0,1
      DO 43 I = 2,N1-L+1
   43 F(I)=Z*F(I-1)+F(I)
      M=M+1
      IF(F(N1) .NE. 0) THEN
       Z1=(ABS(F(N1))*F(N1+1))/F(N1)
       IF(TST(Z1,ZV) .LT. 0) GO TO 44
      ENDIF
   42 CONTINUE
      M1=2
      KK=2*KK
      IF(KK .LE. 8) GO TO 40

   44 S2=1
      S=ABS(F(N1+1))
      S1=ABS(F(N1))
      IF(S .GT. S1*(1+ABS(Z))) S2=S1/(2*S)
      Z=Z-S2*F(N1+1)/F(N1)
      IF(M .GE. MAXIT) THEN
       CALL MTLPRT(NAME,'C209.2','TOO MANY ITERATIONS')
       RETURN
      ENDIF
      DO 81 I = 1,N1+1
   81 F(I)=B(I)
      DO 83 L = 0,2
      DO 83 I = 2,N1-L+1
   83 F(I)=Z*F(I-1)+F(I)
      M=M+1
      GO TO 8

   18 Z=0
      H(K)=0

   13 IF(Z .EQ. 0) GO TO 53
      X=ABS(B(N1+1))
      Z1=1
      JR=0
      DO 51 J = 2,N1
      Z1=Z*Z1
      Y=ABS(Z1*B(N1-J+2))
      IF(Y .GT. X) THEN
       X=Y
       JR=J-1
      ENDIF
   51 CONTINUE
      IF(JR .GT. 0) THEN
       B(N1+1)=-B(N1+1)/Z
       DO 52 J = N1,N1-JR+2,-1
   52  B(J)=(B(J+1)-B(J))/Z
       DO 55 J = N1-JR+1,N1
   55  B(J)=B(J+1)
      ENDIF
      DO 54 J = 2,N1-JR
   54 B(J)=B(J)+Z*B(J-1)

   53 N1=N1-1
      Z0(K)=Z
   10 CONTINUE
      Z0(N)=(-B(2)+SQRT(B(2)**2-4*B(1)*B(3)))/(2*B(1))
      Z0(N-1)=-Z0(N)-B(2)/B(1)
      IF(B(3) .EQ. 0) H(N)=0
      IF(B(2) .EQ. 0 .AND. B(3) .EQ. 0) H(N-1)=0
   12 CONTINUE

      DO 14 K=1,N
      IF(H(K) .EQ. 0) GO TO 14
      H(K)=0

      DO 91 I = 1,N+1
   91 F(I)=C(I-1)
      ZA=ABS(Z0(K))
      HH=HF*ABS(C(0))
      DO 92 I = 2,N+1
      F(I)=Z0(K)*F(I-1)+F(I)
      D=ABS(F(I))
   92 HH=ZA*HH+D
      DO 93 L = 1,2
      DO 93 I = 2,N-L+1
   93 F(I)=Z0(K)*F(I-1)+F(I)
      IF(F(N+1) .EQ. 0 .AND. F(N) .EQ. 0) GO TO 14
      ERR=EPS*(FAC2*HH-D)
      IF(D .LE. ERR) THEN
       S=WN*ERR
      ELSE
       S=WN*ABS(F(N+1))
      ENDIF
      S1=SQRT(ABS(F(N)**2-2*F(N-1)*F(N+1)))
      IF(S .LE. 10*S1*ABS(Z0(K))) THEN
       H(K)=S/S1
      ELSE
       WRITE(ERRTXT,103) M
       CALL MTLPRT(NAME,'C209.3',ERRTXT)
      ENDIF
   14 CONTINUE
      RETURN
103   FORMAT('BOUND FOR Z0',I3,' CANNOT BE FOUND')
      END
