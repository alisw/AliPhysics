*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RCHEBN(M,N,A,MDIM,B,TOL,RELERR,X,RESMAX,IRK,ITER,IOCD)
 
      DIMENSION A(MDIM,*),B(*),X(*)
 
      PARAMETER (R0 = 0, R1 = 1, R2 = 2)
 
      DATA BIG /1E+37/
 
      IRK=N
      NR=1
      RELTMP=RELERR
      RELERR=0
 
      CALL RVSET(M,R1,A(1,N+1),A(2,N+1))
      CALL RVSCL(M,-R1,B(1),B(2),A(1,N+2),A(2,N+2))
      DO 10 I = 1,M
      A(I,N+3)=N+I
   10 CONTINUE
      DO 20 J = 1,N
      A(M+1,J)=J
   20 CONTINUE
      A(M+1,N+1)=0
      ITER=0
      IOCD=1
      CALL RVSET(N,R0,X(1),X(2))
 
      LEV=1
 
      K=0
      NK=N+1
   30 K=K+1
      NK=NK-1
      MODE=0
      CALL RVSET(M-K+1,R1,B(1),B(2))
 
   50 D=-BIG
      DO 60 I = K,M
      IF(B(I) .NE. 0) THEN
       DD=ABS(A(I,N+2))
       IF(DD .GT. D) THEN
        IQ=I
        D=DD
       ENDIF
      ENDIF
   60 CONTINUE
      IF(K .LE. 1 .AND. D .LE. TOL) THEN
       RESMAX=0
       MODE=2
       GOTO 380
      ENDIF
      D=TOL
      DO 80 J = 1,NK
      DD=ABS(A(IQ,J))
      IF(DD .GT. D) THEN
       IP=J
       D=DD
      ENDIF
   80 CONTINUE
      IF(D .GT. TOL) GOTO 330
 
      B(IQ)=0
      IF(MODE .EQ. 1)  GOTO 50
      DO 100 I = K,M
      IF(B(I) .NE. 0) THEN
       DO 90 J = 1,NK
       IF(ABS(A(I,J)) .GT. TOL) THEN
        MODE=1
        GOTO 50
       ENDIF
   90  CONTINUE
      ENDIF
  100 CONTINUE
      IRK=K-1
      NR=N+1-IRK
      IOCD=0
      GOTO 160
 
  110 CALL RVXCH(N+3,A(IQ,1),A(IQ,2),A(K,1),A(K,2))
      CALL RVXCH(M+1,A(1,IP),A(2,IP),A(1,NK),A(2,NK))
      IF(K .LT. N)  GOTO 30
  160 IF(IRK .EQ. M) GOTO 380
 
      LEV=2
 
      D=TOL
      DO 170 I = IRK+1,M
      DD=ABS(A(I,N+2))
      IF(DD .GT. D) THEN
       IQ=I
       D=DD
      ENDIF
  170 CONTINUE
 
      IF(D .LE. TOL) THEN
       RESMAX=0
       MODE=3
       GOTO 380
      ENDIF
      IF(A(IQ,N+2) .GE. -TOL) THEN
       A(IQ,N+1)=2-A(IQ,N+1)
       CALL RVSCL(N+4-NR,-R1,A(IQ,1),A(IQ,2),A(IQ,1),A(IQ,2))
       A(IQ,N+1)=-A(IQ,N+1)
      ENDIF
      DO 220 J = NR,N
      IF(A(IQ,J) .GE. TOL) THEN
       CALL RVSCA(M,R2,
     1            A(1,J),A(2,J),A(1,N+1),A(2,N+1),A(1,N+1),A(2,N+1))
       CALL RVSCL(M,-R1,A(1,J),A(2,J),A(1,J),A(2,J))
       A(M+1,J)=-A(M+1,J)
      ENDIF
  220 CONTINUE
      IP=N+1
      GOTO 330
 
  230 IF(IRK+1 .EQ. M) GO TO 380
      CALL RVXCH(IRK+3,A(IQ,1),A(IQ,2),A(M,1),A(M,2))
 
      LEV=3
 
  260 D=-TOL
      H=2*A(M,N+2)
      DO 280 I = IRK+1,M-1
      IF(A(I,N+2) .LT. D) THEN
       IQ=I
       D=A(I,N+2)
       MODE=0
      ELSE
       DD=H-A(I,N+2)
       IF(DD .LT. D) THEN
        IQ=I
        D=DD
        MODE=1
       ENDIF
      ENDIF
  280 CONTINUE
      IF(D .GE. -TOL)  GOTO 380
      DD=-D/A(M,N+2)
      IF(DD .LT. RELTMP) THEN
       RELERR=DD
       MODE=4
       GOTO 380
      ENDIF
      IF(MODE .NE. 0) THEN
       CALL RVSCS(IRK+1,R2,
     1            A(M,1),A(M,2),A(IQ,1),A(IQ,2),A(IQ,1),A(IQ,2))
       A(IQ,N+2)=D
       A(IQ,N+3)=-A(IQ,N+3)
      ENDIF
      D=BIG
      DO 320 J = NR,N+1
      IF(A(IQ,J) .GT. TOL) THEN
       DD=A(M,J)/A(IQ,J)
       IF(DD .LT. D) THEN
        IP=J
        D=DD
       ENDIF
      ENDIF
  320 CONTINUE
      IF(D .LT. BIG) GO TO 330
      IOCD=2
      GOTO 380
 
  330 RPVT=1/A(IQ,IP)
      CALL RVSCL(M,RPVT,A(1,IP),A(2,IP),A(1,IP),A(2,IP))
      DO 360 I = 1,M
      IF(I .NE. IQ) THEN
       D=A(I,IP)
       CALL RVSCA(N+3-NR,-D,
     1          A(IQ,NR),A(IQ,NR+1),A(I,NR),A(I,NR+1),A(I,NR),A(I,NR+1))
       A(I,IP)=D
      ENDIF
  360 CONTINUE
      CALL RVSCL(IRK+2,-RPVT,A(IQ,NR),A(IQ,NR+1),A(IQ,NR),A(IQ,NR+1))
      A(IQ,IP)=RPVT
      D=A(M+1,IP)
      A(M+1,IP)=A(IQ,N+3)
      A(IQ,N+3)=D
      ITER=ITER+1
      GOTO (110,230,260), LEV
 
  380 CALL RVSET(M,R0,B(1),B(2))
      IF(MODE .EQ. 2)  GOTO 450
      DO 400 I = 1,IRK
      X(INT(A(I,N+3)))=A(I,N+2)
  400 CONTINUE
      IF(MODE .EQ. 3 .OR. IRK .EQ. M)  GOTO 450
      DO 410 J = NR,N+1
      B(INT(ABS(A(M+1,J)))-N)=A(M,N+2)*SIGN(R1,A(M+1,J))
  410 CONTINUE
      DO 420 I = IRK+1,M-1
      B(INT(ABS(A(I,N+3)))-N)=(A(M,N+2)-A(I,N+2))*SIGN(R1,A(I,N+3))
  420 CONTINUE
  430 DO 440 J = NR,N+1
      IF(ABS(A(M,J)) .LE. TOL) THEN
       IOCD=0
       GOTO 450
      ENDIF
  440 CONTINUE
  450 IF(MODE .NE. 2 .AND. MODE .NE. 3) RESMAX=A(M,N+2)
      IF(IRK .EQ. M) RESMAX=0
      IF(MODE .EQ. 4) RESMAX=RESMAX-D
      RETURN
      END
#endif
