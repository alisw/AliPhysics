*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:40  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DLHOIN(A,MA,M,N,MAXV,V,NV,NVEC,EPS,IOUT,W,IW)

C     Solving Systems of Homogeneous Linear Inequalities.
C     Based on K.S. Koelbig and F. Schwarz, A Programm for Solving
C     Systems of Homogeneous Linear Inequalities,
C     Computer Phys. Comm. 17 (1979) 375-382

#include "gen/imp64.inc"
      CHARACTER TEXT*(*)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DLHOIN')
      PARAMETER (TEXT = '+++++ CERN F500 DLHOIN ... ')

      DIMENSION A(MA,*),V(NV,*),IW(MA,*),W(MAXV,*)

      IF(MAXV .LT. N) THEN
       CALL MTLPRT(NAME,'F500.1','MAXV TOO SMALL')
       RETURN
      END IF
C
C*****SETS INITIAL VALUES FOR BOOKKEEPING
C
      DO 3 I = 1,M
      IW(I,1)=0
    3 IW(I,2)=I
      DO 8 I = 1,N
      IW(I,2)=0
    8 IW(I,5)=I
      NINC=N
      NVEC=N
C
C*****DETERMINES N BASIS VECTORS OF THE INITIAL POLYHEDRAL CONE
C
      CALL DMCPY(N,N,A(1,1),A(1,2),A(2,1),V(1,1),V(1,2),V(2,1))
      CALL DINV(N,V,NV,W,IFAIL)
      IF(IFAIL .EQ. -1) THEN
       CALL MTLPRT(NAME,'F500.2','MATRIX A(N,N) SINGULAR')
       RETURN
      END IF
      DO 1 I = 1,NVEC
      S=DVMPY(N,V(1,I),V(2,I),V(1,I),V(2,I))
    1 CALL DVSCL(N,1/SQRT(S),V(1,I),V(2,I),V(1,I),V(2,I))

      IF(IOUT .EQ. 1) THEN
       WRITE(6,111) TEXT
       DO 80 I = 1,NVEC,7
       DO 81 J = 1,N
   81  WRITE(6,'(1X,I9,7E15.6)') J,(V(J,I1),I1=I,MIN(NVEC,I+6))
   80  WRITE(6,'(1X)')
      END IF
C
C*****COMPUTES MATRIX OF SCALAR PRODUCTS
C
   17 DO 20 I = 1,NINC
      DO 20 K = 1,NVEC
   20 W(K,I)=DVMPY(N,A(IW(I,5),1),A(IW(I,5),2),V(1,K),V(2,K))

      IF(IOUT .EQ. 1) THEN
       WRITE(6,112) TEXT
       DO 82 J = 1,NVEC,7
       DO 83 I = 1,NINC
   83  WRITE(6,'(1X,I9,7E15.6)') I,(W(J2,I),J2=J,MIN(NVEC,J+6))
   82  WRITE(6,'(1X)')
      END IF

C*****DETERMINES REDUNDANT INEQUALITIES AND CHOOSES NEW ONE

      DO 40 K = 1,M
      IW(K,3)=0
      IF(IW(K,2) .EQ. 0) GO TO 40
      DO 45 I = 1,NVEC
      IF(DVMPY(N,A(K,1),A(K,2),V(1,I),V(2,I)) .GT. 0) IW(K,3)=IW(K,3)+1
   45 CONTINUE
   40 CONTINUE
      NNEG=NVEC
      DO 48 K = 1,M
      IF(IW(K,2) .EQ. 0) GO TO 48
      IF(IW(K,3) .EQ. 0) THEN
       WRITE(ERRTXT,103) K
       CALL MTLPRT(NAME,'F500.3',ERRTXT)
       RETURN
      END IF
      IF(IW(K,3) .EQ. NVEC) THEN
       IW(K,1)=K
       IW(K,2)=0
      END IF
      IF(IW(K,3) .LT. NNEG) THEN
       NNEG=IW(K,3)
       KNEW=K
      END IF
   48 CONTINUE
      IF(NNEG .EQ. NVEC) THEN
       DO 74 I = 1,M
       IW(I,1)=I
       DO 75 J = 1,NVEC
       IF(DVMPY(N,A(I,1),A(I,2),V(1,J),V(2,J)) .GE. EPS) GO TO 75
       IW(I,1)=0
       GO TO 74
   75  CONTINUE
   74  CONTINUE
       RETURN
      END IF

      IF(IOUT .EQ. 1) WRITE(6,113) TEXT,KNEW
      IW(KNEW,2)=0
C
C*****COMPUTES VECTOR OF SCALAR PRODUCTS
C
      DO 50 I = 1,NVEC
   50 W(I,M+1)=DVMPY(N,A(KNEW,1),A(KNEW,2),V(1,I),V(2,I))

      IF(IOUT .EQ. 1) THEN
       WRITE(6,114) TEXT
       DO 84 J = 1,NVEC,7
   84  WRITE(6,'(1X,I9,7E15.6)') J,(W(J2,M+1),J2=J,MIN(NVEC,J+6))
       WRITE(6,'(1X)')
      END IF
C
C*****DETERMINES BASIS VECTORS FOR NEW CONE
C
      NTVE=NVEC
      DO 60 I = 1,NVEC-1
      DO 60 J = I+1,NVEC
      IF(W(I,M+1)*W(J,M+1) .GT. 0) GO TO 60
      NT=0
      DO 62 L = 1,NINC
      IW(L,4)=1
      IF(ABS(W(I,L)) .GT. EPS .OR. ABS(W(J,L)) .GT. EPS) GO TO 62
      NT=NT+1
      IW(L,4)=0
   62 CONTINUE
      IF(NT .LT. N-2) GO TO 60
      DO 63 K = 1,NVEC
      IF(K .EQ. I .OR. K .EQ. J) GO TO 63
      MT=0
      DO 64 L = 1,NINC
      IF(IW(L,4) .EQ. 0 .AND. ABS(W(K,L)) .LT. EPS) MT=MT+1
   64 CONTINUE
      IF(MT .EQ. N-2) GO TO 60
   63 CONTINUE
      NTVE=NTVE+1
      IF(NTVE .GT. MAXV) THEN
       CALL MTLPRT(NAME,'F500.1','MAXV TOO SMALL')
       RETURN
      END IF
      DO 65 L = 1,N
   65 V(L,NTVE)=ABS(W(J,M+1))*V(L,I)+ABS(W(I,M+1))*V(L,J)
   60 CONTINUE
      DO 66 I = 1,NTVE
      S=DVMPY(N,V(1,I),V(2,I),V(1,I),V(2,I))
   66 CALL DVSCL(N,1/SQRT(S),V(1,I),V(2,I),V(1,I),V(2,I))
C
C*****ELIMINATES VECTORS WITH NEGATIVE SCALAR PRODUCT
C
      NNEW=0
      DO 70 I = 1,NVEC
      IF(W(I,M+1) .LT. 0) GO TO 70
      NNEW=NNEW+1
      CALL DVCPY(N,V(1,I),V(2,I),V(1,NNEW),V(2,NNEW))
   70 CONTINUE
      DO 71 I = NVEC+1,NTVE
      NNEW=NNEW+1
   71 CALL DVCPY(N,V(1,I),V(2,I),V(1,NNEW),V(2,NNEW))
      CALL DMSET(N,NTVE-NNEW,0D0,V(1,NNEW+1),V(1,NNEW+2),V(2,NNEW+1))
      NVEC=NNEW
      NINC=NINC+1
      IW(NINC,5)=KNEW

      IF(IOUT .EQ. 1) THEN
       WRITE(6,115) TEXT
       DO 86 I = 1,NVEC,7
       DO 87 J = 1,N
   87  WRITE(6,'(1X,I9,7E15.6)') J,(V(J,I2),I2=I,MIN(NVEC,I+6))
   86  WRITE(6,'(1X)')
      END IF
      GO TO 17
  103 FORMAT('INEQUALITY ',I5,' IS INCONSISTENT')
  111 FORMAT(7X,A27,'THE N BASIS VECTORS OF THE FIRST CONE'/)
  112 FORMAT(7X,A27,'THE MATRIX OF SCALAR PRODUCTS'/)
  113 FORMAT(7X,A27,'THE NEW INEQUALITY HAS INDEX',I5/)
  114 FORMAT(7X,A27,'THE VECTOR OF SCALAR PRODUCTS'/)
  115 FORMAT(7X,A27,'THE MATRIX OF BASIS VECTORS'/)
      END
#endif
