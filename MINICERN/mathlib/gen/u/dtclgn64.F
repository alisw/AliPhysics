*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DTCLGN(J1,J2,J3,M1,M2,M3,DNUM,DDEN,KPEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DTCLGN(J1,J2,J3,M1,M2,M3,DNUM,DDEN,KPEX)
#include "gen/imp64.inc"
#endif
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DTCLGN')
 
      PARAMETER (MXP = 40)
 
      DIMENSION KPEX(*)
      DIMENSION IPRIM(MXP),PRIM(MXP),LS(MXP,12),LR(MXP,100)
      DIMENSION ISD(11),JS(100),IA(7)
 
      DATA (IA(I),I=1,7) /-1,1,1,-1,-1,1,1/
 
      DATA (IPRIM(I),I=1,MXP)
     1/  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
     2  59, 61, 67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,127,131,
     3 137,139,149,151,157,163,167,173/
 
      DO 11 J = 1,MXP
      LS(J,12)=0
      PRIM(J)=IPRIM(J)
   11 KPEX(J)=0
      DNUM=0
      DDEN=1
      IF(J1 .LT. 0 .OR. J2 .LT. 0 .OR. J3 .LT. 0) GO TO 99
      IF(MOD(J1+J2+J3,2) .NE. 0) GO TO 99
      MX=ABS(M1)
      MY=ABS(M2)
      MZ=ABS(M3)
      IF(J1 .LT. MX .OR. J2 .LT. MY .OR. J3 .LT. MZ) GO TO 99
      IF(MOD(J1+MX,2) .NE. 0 .OR. MOD(J2+MY,2) .NE. 0) GOTO 99
      IF(MOD(J3+MZ,2) .NE. 0) GOTO 99
      IF(M1+M2 .NE. M3) GO TO 99
      ISD(1)=J3+1
      ISD(2)=(J3+J1-J2)/2
      ISD(3)=(J3-J1+J2)/2
      ISD(4)=(J1+J2-J3)/2
      ISD(5)=(J3-M3)/2
      ISD(6)=(J3+M3)/2
      ISD(7)=(J1+J2+J3)/2+1
      ISD(8)=(J1-M1)/2
      ISD(9)=(J1+M1)/2
      ISD(10)=(J2-M2)/2
      ISD(11)=(J2+M2)/2
 
      IF(ISD(2) .LT. 0 .OR. ISD(3) .LT. 0 .OR. ISD(4) .LT. 0) GO TO 99
 
      DO 41 K = 1,MXP
      DO 41 L = 1,11
   41 LS(K,L)=0
      IF(ISD(1) .GT. 1) THEN
       KOVF=1
       N1=ISD(1)
       K=0
   42  K=K+1
       IF(K .EQ. MXP+1) GO TO 48
   43  IF(MOD(N1,IPRIM(K)) .NE. 0) GO TO 42
       LS(K,1)=LS(K,1)+1
       N1=N1/IPRIM(K)
       IF(N1 .NE. 1) GO TO 43
   48  KOVF=MAX(KOVF,N1)
       IF(KOVF .NE. 1) THEN
        WRITE(ERRTXT,101) J1,J2,J3,M1,M2,M3
        CALL MTLPRT(NAME,'U112.1',ERRTXT)
        GO TO 99
       ENDIF
      ENDIF
      DO 59 L = 2,11
      DO 54 N = 2,ISD(L)
      KOVF=1
      N1=N
      K=0
   52 K=K+1
      IF(K .EQ. MXP+1) GO TO 58
   53 IF(MOD(N1,IPRIM(K)) .NE. 0) GO TO 52
      LS(K,L)=LS(K,L)+1
      N1=N1/IPRIM(K)
      IF(N1 .NE. 1) GO TO 53
   58 KOVF=MAX(KOVF,N1)
      IF(KOVF .NE. 1) THEN
       WRITE(ERRTXT,101) J1,J2,J3,M1,M2,M3
       CALL MTLPRT(NAME,'U112.1',ERRTXT)
       GO TO 99
      ENDIF
   54 CONTINUE
   59 CONTINUE
 
      DO 5 J = 1,MXP
    5 LS(J,12)=LS(J,12)+LS(J,1)+LS(J,2)+LS(J,3)+LS(J,4)+LS(J,5)+LS(J,6)
     1        -LS(J,7)-LS(J,8)-LS(J,9)-LS(J,10)-LS(J,11)
      MM=0
      ISD(1)=J2+J3+M1
      ISD(2)=J1-M1
      ISD(3)=0
      ISD(4)=J3-J1+J2
      ISD(5)=J3+M3
      ISD(6)=J1-J2-M3
      ISD(7)=J2+M2
      DO 12 L = 1,7
   12 ISD(L)=ISD(L)/2-IA(L)
    3 DO 13 L = 1,7
   13 ISD(L)=ISD(L)+IA(L)
      IF(ISD(2) .LT. 0 .OR. ISD(6) .LT. 0) GO TO 3
      IF(ISD(1) .LT. 0 .OR. ISD(4) .LT. 0 .OR. ISD(5) .LT. 0) GO TO 4
 
      DO 69 L = 1,6
      DO 61 K = 1,MXP
   61 LS(K,L)=0
      DO 64 N = 2,ISD(L)
      KOVF=1
      N1=N
      K=0
   62 K=K+1
      IF(K .EQ. MXP+1) GO TO 68
   63 IF(MOD(N1,IPRIM(K)) .NE. 0) GO TO 62
      LS(K,L)=LS(K,L)+1
      N1=N1/IPRIM(K)
      IF(N1 .NE. 1) GO TO 63
   68 KOVF=MAX(KOVF,N1)
      IF(KOVF .NE. 1) THEN
       WRITE(ERRTXT,101) J1,J2,J3,M1,M2,M3
       CALL MTLPRT(NAME,'U112.1',ERRTXT)
       GO TO 99
      ENDIF
   64 CONTINUE
   69 CONTINUE
 
      MM=MM+1
      JS(MM)=(-1)**ISD(7)
      DO 9 J = 1,MXP
    9 LR(J,MM)=LS(J,1)+LS(J,2)-LS(J,3)-LS(J,4)-LS(J,5)-LS(J,6)
      GO TO 3
 
    4 DO 31 J = 1,MXP
      LS(J,8)=10000
      DO 31 JK = 1,MM
   31 LS(J,8)=MIN(LS(J,8),LR(J,JK))
      BSUM=0
      DO 32 JM = 1,MM
      BNUM=1
      DO 33 J = 1,MXP
      JEX=LR(J,JM)-LS(J,8)
      IF(JEX .GT. 0) BNUM=BNUM*PRIM(J)**JEX
   33 CONTINUE
   32 BSUM=BSUM+BNUM*JS(JM)
      IF(BSUM .EQ. 0) GO TO 99
      DO 71 K = 1,MXP
   71 LS(K,1)=0
      ASUM=ABS(BSUM)
      IF(ASUM .GT. 1) THEN
       OVF=1
       A1=ASUM
       K=0
   72  K=K+1
       IF(K .EQ. MXP+1) GO TO 78
   73  IF(MOD(A1,PRIM(K)) .NE. 0) GO TO 72
       LS(K,1)=LS(K,1)+1
       A1=A1/PRIM(K)
       IF(A1 .NE. 1) GO TO 73
   78  OVF=MAX(OVF,A1)
       IF(OVF .NE. 1) THEN
        WRITE(ERRTXT,101) J1,J2,J3,M1,M2,M3
        CALL MTLPRT(NAME,'U112.1',ERRTXT)
        GO TO 99
       ENDIF
      ENDIF
      DO 22 J = 1,MXP
   22 KPEX(J)=LS(J,12)+2*(LS(J,8)+LS(J,1))
      BNUM=1
      BDEN=1
      DO 16 J = 1,MXP
      PP=PRIM(J)**ABS(KPEX(J))
      IF(KPEX(J) .GE. 0) THEN
       BNUM=BNUM*PP
      ELSE
       BDEN=BDEN*PP
      ENDIF
   16 CONTINUE
      DNUM=SIGN(BNUM,BSUM)
      DDEN=BDEN
   99 RETURN
  101 FORMAT('LIST OF PRIME NUMBERS EXHAUSTED',2X,3I5,2X,3I5)
      END
