*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MINSQ (M,N,F,X,E,IPRINT,NFUN ,NW,W,COV,XSTEP)
      DIMENSION F(M),X(N),E(N),W(NW)
      IFLAG=1
      ESCALE=1.
      MAXFUN=NFUN
      MPLUSN=M+N
      KST=N+MPLUSN
      NPLUS=N+1
      KINV=NPLUS*(MPLUSN+1)
      KSTORE=KINV-MPLUSN-1
      CALL FCN1(M,N,F,X,IFLAG,NW,W,KEND,KEND2N)
      IFLAG=4
      NN=N+N
      K=NN
      DO 1 I=1,M
      K=K+1
      W(K)=F(I)
    1 CONTINUE
      IINV=2
      K=KST
      I=1
    2 X(I)=X(I)+E(I)
      CALL FCN1(M,N,F,X,IFLAG,NW,W,KEND,KEND2N)
      X(I)=X(I)-E(I)
      DO 3 J=1,N
      K=K+1
      W(K)=0.
      W(J)=0.
    3 CONTINUE
      SUM=0.
      KK=NN
      DO 4 J=1,M
      KK=KK+1
      F(J)=F(J)-W(KK)
      SUM=SUM+F(J)*F(J)
    4 CONTINUE
      IF (SUM) 5,5,6
    5 WRITE(6,7)I
    7 FORMAT(5X,'E(',I3,') UNREASONABLY SMALL')
      RETURN
    6 SUM=1./SQRT (SUM)
      J=K-N+I
      W(J)=E(I)*SUM
      DO 9 J=1,M
      K=K+1
      W(K)=F(J)*SUM
      KK=NN+J
      DO 11 II=1,I
      KK=KK+MPLUSN
      W(II)=W(II)+W(KK)*W(K)
   11 CONTINUE
    9 CONTINUE
      ILESS=I-1
      IGAMAX=N+I-1
      INCINV=N-ILESS
      INCINP=INCINV+1
      IF (ILESS) 13,13,14
   13 W(KINV)=1.
      GO TO 15
   14 B=1.
      DO 16 J=NPLUS,IGAMAX
      W(J)=0.
   16 CONTINUE
      KK=KINV
      DO 17 II=1,ILESS
      IIP=II+N
      W(IIP)=W(IIP)+W(KK)*W(II)
      JL=II+1
      IF (JL-ILESS) 18,18,19
   18 DO 20 JJ=JL,ILESS
      KK=KK+1
      JJP=JJ+N
      W(IIP)=W(IIP)+W(KK)*W(JJ)
      W(JJP)=W(JJP)+W(KK)*W(II)
   20 CONTINUE
   19 B=B-W(II)*W(IIP)
      KK=KK+INCINP
   17 CONTINUE
      B=1./B
      KK=KINV
      DO 21 II=NPLUS,IGAMAX
      BB=-B*W(II)
      DO 22 JJ=II,IGAMAX
      W(KK)=W(KK)-BB*W(JJ)
      KK=KK+1
   22 CONTINUE
      W(KK)=BB
      KK=KK+INCINV
   21 CONTINUE
      W(KK)=B
   15 GO TO (27,24),IINV
   24 I=I+1
      IF (I-N) 2,2,25
   25 IINV=1
      FF=0.
      KL=NN
      DO 26 I=1,M
      KL=KL+1
      F(I)=W(KL)
      FF=FF+F(I)*F(I)
   26 CONTINUE
      ICONT = 1
      ISS=1
      MC=N+1
      IPP=IPRINT*(IPRINT-1)
      ITC=0
      IPS=1
      IPC=0
   27 IPC=IPC-IPRINT
      IF (IPC) 28,29,29
   28 WRITE(6,30)ITC,MC,FF
   30 FORMAT(//5X,'ITERATION',I4,I9,' CALLS OF FCN   ',5X,'F=',E24.14)
      WRITE(6,31)(X(I),I=1,N)
   31 FORMAT(5X,'VARIABLES',/(5E24.14))
      WRITE(6,32)(F(I),I=1,M)
   32 FORMAT(5X,'FUNCTIONS',/(5E24.14))
      IPC=IPP
      GO TO (29,33),IPS
   29 GO TO (34,35),ICONT
   35 IF (CHANGE-1.) 10,10,36
   10 CONTINUE
C*UL   37 WRITE(6,38)
      WRITE(6,38)
   38 FORMAT(//5X,'      FINAL VALUES OF FUNCTIONS AND VARIABLES')
      IPS=2
      GO TO 28
   33 IFLAG=3
      CALL FCN1(M,N,F,X,IFLAG,NW,W,KEND,KEND2N)
      IF(COV.EQ.0.) RETURN
      IWC1 = 2*N+N*(N+M)
      IWC=IWC1
      DO 91 I=1,N
      DO 91 J=1,I
      IWC=IWC+1
      W(IWC)=0.
      DO 90 MR=1,N
      JCOV=N+J+MR*(N+M)
      ICOV=N+I+MR*(N+M)
   90 W(IWC)=W(IWC)+W(ICOV)*W(JCOV)
   91 CONTINUE
      WRITE(6,1000)
 1000 FORMAT(///40X,'VARIANCE-COVARIANCE MATRIX'///)
      IB=IWC1
      DO 92 I=1,N
      IA=IB+1
      IB=IB+I
   92 WRITE(6,1001)(W(J),J=IA,IB)
 1001 FORMAT(10X,5G20.6//)
      RETURN
   36 ICONT=1
   34 ITC=ITC+1
      K=N
      KK=KST
      DO 39 I=1,N
      K=K+1
      W(K)=0.
      KK=KK+N
      W(I)=0.
      DO 40 J=1,M
      KK=KK+1
      W(I)=W(I)+W(KK)*F(J)
   40 CONTINUE
   39 CONTINUE
      DM=0.
      K=KINV
      DO 41 II=1,N
      IIP=II+N
      W(IIP)=W(IIP)+W(K)*W(II)
      JL=II+1
      IF (JL-N) 42,42,43
   42 DO 44 JJ=JL,N
      JJP=JJ+N
      K=K+1
      W(IIP)=W(IIP)+W(K)*W(JJ)
      W(JJP)=W(JJP)+W(K)*W(II)
   44 CONTINUE
      K=K+1
   43 IF (DM-ABS (W(II)*W(IIP))) 45,41,41
   45 DM=ABS (W(II)*W(IIP))
      KL=II
   41 CONTINUE
      II=N+MPLUSN*KL
      CHANGE=0.
      DO 46 I=1,N
      JL=N+I
      W(I)=0.
      DO 47 J=NPLUS,NN
      JL=JL+MPLUSN
      W(I)=W(I)+W(J)*W(JL)
   47 CONTINUE
      II=II+1
      W(II)=W(JL)
      W(JL)=X(I)
      IF (ABS (E(I)*CHANGE)-ABS (W(I))) 48,48,46
   48 CHANGE=ABS (W(I)/E(I))
   46 CONTINUE
      DO 49 I=1,M
      II=II+1
      JL=JL+1
      W(II)=W(JL)
      W(JL)=F(I)
   49 CONTINUE
      FC=FF
      ACC=0.1/CHANGE
      IT=3
      XC=0.
      XL=0.
      IS=3
      IF (CHANGE-1.) 50,50,51
   50 ICONT=2
   51 CALL VD01A (IT,XC,FC,20,ACC,0.1,XSTEP)
      GO TO (52,53,53,53),IT
   52 MC=MC+1
      IF (MC-MAXFUN) 54,54,55
   55 WRITE(6,56)MAXFUN
   56 FORMAT(5X,I6,' CALLS OF FCN')
      ISS=2
      GO TO 53
   54 XL=XC-XL
      DO 57 J=1,N
      X(J)=X(J)+XL*W(J)
   57 CONTINUE
      XL=XC
      CALL FCN1(M,N,F,X,IFLAG,NW,W,KEND,KEND2N)
      FC=0.
      DO 58 J=1,M
      FC=FC+F(J)*F(J)
   58 CONTINUE
      GO TO (59,59,60),IS
   60 K=N
      IF (FC-FF) 61,51,62
   61 IS=2
      FMIN=FC
      FSEC=FF
      GO TO 63
   62 IS=1
      FMIN=FF
      FSEC=FC
      GO TO 63
   59 IF (FC-FSEC) 64,51,51
   64 K=KSTORE
      GO TO (75,74),IS
   75 K=N
   74 IF (FC-FMIN) 65,51,66
   66 FSEC=FC
      GO TO 63
   65 IS=3-IS
      FSEC=FMIN
      FMIN=FC
   63 DO 67 J=1,N
      K=K+1
      W(K)=X(J)
   67 CONTINUE
      DO 68 J=1,M
      K=K+1
      W(K)=F(J)
   68 CONTINUE
      GO TO 51
   53 K=KSTORE
      KK=N
      GO TO (69,70,69),IS
   70 K=N
      KK=KSTORE
   69 SUM=0.
      DM=0.
      JJ=KSTORE
      DO 71 J=1,N
      K=K+1
      KK=KK+1
      JJ=JJ+1
      X(J)=W(K)
      W(JJ)=W(K)-W(KK)
   71 CONTINUE
      DO 72 J=1,M
      K=K+1
      KK=KK+1
      JJ=JJ+1
      F(J)=W(K)
      W(JJ)=W(K)-W(KK)
      SUM=SUM+W(JJ)*W(JJ)
      DM=DM+F(J)*W(JJ)
   72 CONTINUE
      GO TO (73,10),ISS
   73 J=KINV
      KK=NPLUS-KL
      DO 76 I=1,KL
      K=J+KL-I
      J=K+KK
      W(I)=W(K)
      W(K)=W(J-1)
   76 CONTINUE
      IF (KL-N) 77,78,78
   77 KL=KL+1
      JJ=K
      DO 79 I=KL,N
      K=K+1
      J=J+NPLUS-I
      W(I)=W(K)
      W(K)=W(J-1)
   79 CONTINUE
      W(JJ)=W(K)
      B=1./W(KL-1)
      W(KL-1)=W(N)
      GO TO 88
   78 B=1./W(N)
   88 K=KINV
      DO 80 I=1,ILESS
      BB=B*W(I)
      DO 81 J=I,ILESS
      W(K)=W(K)-BB*W(J)
      K=K+1
   81 CONTINUE
      K=K+1
   80 CONTINUE
      IF (FMIN-FF) 82,83,83
   83 CHANGE=0.
      GO TO 84
   82 FF=FMIN
      CHANGE=ABS (XC)*CHANGE
   84 XL=-DM/FMIN
      SUM=1./SQRT (SUM+DM*XL)
      K=KSTORE
      DO 85 I=1,N
      K=K+1
      W(K)=SUM*W(K)
      W(I)=0.
   85 CONTINUE
      DO 86 I=1,M
      K=K+1
      W(K)=SUM*(W(K)+XL*F(I))
      KK=NN+I
      DO 87 J=1,N
      KK=KK+MPLUSN
      W(J)=W(J)+W(KK)*W(K)
   87 CONTINUE
   86 CONTINUE
      GO TO 14
      END
