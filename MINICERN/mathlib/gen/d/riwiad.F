*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:16  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RIWIAD(F)
#include "gen/imp64.inc"
C
C********** COMMUNICATION TO CALLING PROGRAM
C
      COMMON/PARAMS/ACC,NDIM,NSUB,ITER
      COMMON/ANSWER/VALUE,ERROR
      COMMON/STORE/XA(11),XB(11),XC(11),XD(11),MA(11),MB(11),MC(11)
      COMMON/STORE1/R(10000),LR
      COMMON/OPTION/IPRINT,ICONV,IRESET
      COMMON/RANDOM/NSHOTS
      COMMON/INTERN/FACTOR,ALFA,BETA,GAMMA,DELTA,LEVEL,NMIN
C*NS  INTEGER XLOC, FNAME
      EXTERNAL F
C
#if defined(CERNLIB_DOUBLE)
      REAL Q
C     MODIFICATIONS FOR DOUBLE PRECISION
      REAL(N)=DBLE(N)
C............
#endif
      CALL RIWIBD
C***************CONSTANTS AND COUNTERS
      ZERO=1.0E-12
      ONE=1.0-ZERO
      ITC=0
      NCALL=0
      SWI=0.0
      SRV=0.0
      SWISQ=0.0
      ANMIN=NMIN
      ANDIM=NDIM
      ANSUB=NSUB
      RNDIM=1.0/ANDIM
      ANSHT1=NSHOTS-1
      RNSHTS=1.0/NSHOTS
      RNSHT1=1.0/(NSHOTS-1)
C
C********** INITIALISE INTERVALS
C
      IF(MA(1).EQ.0) GOTO 201
      IF(IRESET.EQ.0) GOTO 202
  201 K=INT(ANSUB**RNDIM+ZERO)
      IF(K.LT.NMIN) GOTO 207
      A=1.0/K
      L=NDIM*K
      DO 204 I=1,L
  204 R(I)=A
      DO 203 I=1,NDIM
  203 MA(I)=K
  202 CONTINUE
C
C********** PRINT OUT INITIAL VALUES
C
      IF(IPRINT.EQ.0) GOTO 1
      WRITE(6,2)
    2 FORMAT('1........ RIEMANN INTEGRATION WITH INTERVAL ADJUSTMENT',
     *'........')
      WRITE(6,3)      ACC,NDIM,NSUB,ITER
3      FORMAT('0******** INPUT PARAM ********'//
     +,' NAME OF INTEGRAND...',13X,'F   '/
     +,' RELATIVE ACCURACY... ',12X,F7.5/
     +,' DIMENSION OF INTEGRAL... ',I15/
     +,' MAXIMUM NUMBER OF SUBVOLUMES... ',I8/
     +,' MAXIMUM NUMBER OF ITERATIONS... ',I8)
      IF(IPRINT.EQ.2) GOTO 1
      WRITE(6,4)NSHOTS,LEVEL,FACTOR,IPRINT,
     1ICONV,IRESET,ALFA,BETA,GAMMA,DELTA,NMIN
    4 FORMAT('0'/'0******** INTERNAL PARAMETERS ********'//
     +,' NUMBER OF CALLS PER SUBVOLUME... ',I7/
     +,' LEVEL OF CONFIDENCE (PER CENT)...',I7/
     +,' CORRESPONDING FACTOR...',F17.2/
     +,' PRINT OPTION...',I25/
     +,' CONVERGENCE CRITERION...',I16/
     +,' RESET OPTION...',I25/
     +,' RATE OF CONVERGENCE...',2X,4F4.2/
     +,' MINIMUM NUMBER OF INTERVALS...',I10)
C
C********** ENTRY AND REENTRY FOR ITERATIONS
C
    1 ITC=ITC+1
      J=1
      NID114=0
      DO 6 I=1,NDIM
      MB(I)=NID114
      K=MA(I)
      J=J*K
    6 NID114=NID114+K
      NCALL=NCALL+NSHOTS*J
      LA=0
      LB=NID114
      LC=2*NID114
      LD=3*NID114
      LE=4*NID114
      IF(LE.GT.LR) GOTO 107
C
C********** INITIALISE INTEGRATION LOOP
C
C*UL7 L=0
      L=0
      V=1.0
      DO 110 I=1,NDIM
      MC(I)=1
      K=L+1
      L=L+MA(I)
      A=R(LA+K)
      XA(I)=A
      XB(I)=0.0
      XC(I)=V
      V=V*A
      X=0.0
      DO 8 J=K,L
      R(LB+J)=X
      X=X+R(LA+J)
      R(LC+J)=0.0
    8 R(LD+J)=0.0
      R(LA+L)=R(LA+L)+1.0-X
  110 CONTINUE
C
C********** INTEGRATION LOOP
C
   15 SF1=0.0
      SF2=0.0
      DO 9 N=1,NSHOTS
      DO 10 I=1,NDIM
   10 XD(I)=XB(I)+XA(I)*RNDM(I)
      Y=F(XD)
      SF1=SF1+Y
    9 SF2=SF2+Y**2
      SF1=SF1*RNSHTS
      SF2=SF2*RNSHTS
      C=V*SF1
      D=V**2*(SF2-SF1**2)*RNSHT1
      DO 11 I=1,NDIM
      J=MB(I)+MC(I)
      R(LC+J)=R(LC+J)+C
   11 R(LD+J)=R(LD+J)+D
C
C********** LOOP THE LOOP
C
      I=NDIM
   13 J=MC(I)
      K=MB(I)+1
      IF(J.NE.MA(I)) GOTO 12
      MC(I)=1
      XA(I)=R(LA+K)
      XB(I)=0.0
      I=I-1
      IF(I) 14,14,13
   12 MC(I)=J+1
      V=XC(I)
      K=K+J
      A=R(LA+K)
      XA(I)=A
      XB(I)=R(LB+K)
      V=V*A
   16 IF(I.EQ.NDIM) GOTO 15
      I=I+1
      XC(I)=V
      V=V*XA(I)
      GOTO 16
C
C********** CALCULATE INTEGRAL AND ERROR
C
   14 J=MA(1)
      SI=0.0
      SV=0.0
      DO 17 I=1,J
      SI=SI+R(LC+I)
   17 SV=SV+R(LD+I)
      IF(SV.EQ.0.0) GOTO 400
      RV=1.0/SV
      WI=SI*RV
      SRV=SRV+RV
      SWI=SWI+WI
      VALUE=SWI/SRV
      FSV=FACTOR*SQRT(SV)
      ERROR=FACTOR/SQRT(SRV)
      SWISQ=SWISQ+RV*SI**2
      CHISQ=SWISQ-VALUE**2*SRV
C
C********** REDUCE TO RELATIVE QUANTITIES
C
      RI=1.0/SI
      DO 18 J=1,NID114
      R(LC+J)=R(LC+J)*RI
   18 R(LD+J)=R(LD+J)*RV
C
C********** PRINT OUT THE RESULTS OF THIS ITERATION
C
      IF(IPRINT.EQ.0) GOTO 25
      WRITE(6,20)ITC,SI,FSV,VALUE,ERROR,CHISQ
     1,NCALL,(MA(I),I=1,NDIM)
 20   FORMAT('0'//////' ******** ITERATION NUMBER',I3,' ********'//
     +,' VALUE OF INTEGRAL...',E13.5/
     +,' ESTIMATED ERROR.....',E13.5//
     +,' ACCUMULATED VALUE...',E13.5/
     +,' ACCUMULATED ERROR...',E13.5//
     +,' CHI SQUARE..........',F13.2//
     +,' NUMBER OF CALLS TO F',I13/
     +,' INTERVAL STRUCTURE..',5X,(25I4)//)
      IF(IPRINT.EQ.2) GOTO 25
      DO 19 I=1,NDIM
      WRITE(6,21)I
   21 FORMAT('0'/'0........ AXIS NUMBER',I3)
      K=MB(I)
      L=MA(I)+K
      K=K+1
      WRITE(6,22)(R(LA+J),J=K,L)
   22 FORMAT('0INTERVAL CHOICE'/(10E13.2))
      WRITE(6,23)(R(LC+J),J=K,L)
   23 FORMAT('0RELATIVE CONTRIBUTION TO INTEGRAL'/(10E13.2))
      WRITE(6,24)(R(LD+J),J=K,L)
   24 FORMAT('0RELATIVE CONTRIBUTION TO ERROR'/(10E13.2))
   19 CONTINUE
C
C********** CHECK IF ACCURACY REACHED OR END OF ITERATIONS
C
   25 ACC1=ERROR/ABS(VALUE)
      IF(ACC1.GT.ACC) GOTO 26
   28 IF(IPRINT.NE.0) WRITE(6,27) 'WAS',ITC,VALUE,ERROR
      RETURN
  26  IF(ITC.EQ.ITER) THEN
        IF(IPRINT.NE.0) WRITE(6,27)'NOT',ITC,VALUE,ERROR
        RETURN
      ENDIF
C
C********** RECOMPUTE INTERVAL SIZES
C
      L=0
      V=1.0
      DO 29 I=1,NDIM
      K=L+1
      M=MA(I)
      L=L+M
      X=0.0
      Y=0.0
      Z=0.0
      DO 30 J=K,L
      A=R(LA+J)
      C=R(LC+J)
      D=R(LD+J)
      GOTO (301,302,303,304),ICONV
  301 CONTINUE
      P=(1.0-A+ZERO)**ALFA
      Q=(-LOG(D+ZERO)+ZERO)**BETA
      GOTO 399
  302 CONTINUE
      P=A**ALFA
      Q=(D+ZERO)**(-BETA)
      GOTO 399
  303 CONTINUE
      P=((-LOG(D+ZERO)+ZERO)/(1.0-D+ZERO))**ALFA
      Q=((-LOG(A+ZERO)+ZERO)/(1.0-A+ZERO))**(-BETA)
      GOTO 399
  304 CONTINUE
      P=(1.0-LOG(D+ZERO))**ALFA
      Q=(1.0-LOG(A+ZERO))**(-BETA)
      GOTO 399
  399 CONTINUE
      B=A*P*Q
      R(LB+J)=B
      X=X+B
      Y=Y+D/A
   30 Z=Z+C**2/A
      B=ANSHT1 *Y*SV+(Z-1.0)*SI**2
      Y=B**GAMMA*M**DELTA
      V=V*Y
      XA(I)=Y
      X=1.0/X
      MC(I)=0
      DO 31 J=K,L
   31 R(LB+J)=R(LB+J)*X
   29 CONTINUE
C
C********** RECOMPUTE INTERVAL NUMBERS
C
      L=0
   34 Q=(ANSUB/(ANMIN**L*V))**(1.0/(NDIM-L))
      V=1.0
      K=0
      DO 32 I=1,NDIM
      M=MC(I)
      IF(M.EQ.NMIN) GOTO 32
      A=XA(I)*Q
      M=INT(A+0.5)
      IF(M.LE.NMIN) GOTO 33
      XA(I)=A
      V=V*A
      MC(I)=M
      GOTO 32
   33 MC(I)=NMIN
      K=K+1
   32 CONTINUE
      L=L+K
      IF(K.NE.0) GOTO 34
C
C********** COMBINE SIZES AND NUMBERS INTO NEW INTERVALS
C
      NID114=0
      DO 35 I=1,NDIM
   35 NID114=NID114+MC(I)
      IF(5*NID114.GT.LR) GOTO 107
      L=0
      L1=0
      DO 36 I=1,NDIM
      K=L+1
      K1=L1+1
      M=MC(I)
      M1=MA(I)
      L=L+M
      L1=L1+M1
      A=REAL(M1)/REAL(M)
      A=A*ONE
      C=0.0
      DO 38 J=K,L
      B=C
      C=C+A
      I1=INT(B)
      I2=INT(C)
      X=0.0
      DO 39 I3=I1,I2
   39 X=X+R(LB+K1+I3)
   38 R(LE+J)=X-(I2+1-C)*R(LB+K1+I2)-(B-I1)*R(LB+K1+I1)
   36 CONTINUE
      DO 40 L=1,NID114
   40 R(LA+L)=R(LE+L)
      DO 41 I=1,NDIM
   41 MA(I)=MC(I)
      GOTO 1
C
C********** NOT STORAGE ENOUGH
C
  107 X1=ANDIM*ANSUB**RNDIM
      X2=ANMIN*(ANDIM-1.0)+ANSUB/ANMIN**(ANDIM-1.0)
      X1=5.0*X1
      X2=5.0*X2
      X3=SQRT(X1*X2)
      I1=INT(X1)
      I2=INT(X2)
      I3=INT(X3)
      WRITE(6,108)LR,I1,I2,I3
  108 FORMAT('0'/' STORAGE LENGTH =',I5,'  IS TOO SMALL'/
     +,' MINIMUM NECESSARY -',I7/' MAXIMUM POSSIBLE - ',I7/
     +,' PROBABLE VALUE - ',I9)
      STOP 1
C
C********** TOO SMALL NUMBER OF SUBVOLUMES
C
  207 L=NMIN**NDIM
      WRITE(6,206)NSUB,L
  206 FORMAT('0 NSUB TOO SMALL...',I7/
     +,' MINIMUM NECESSARY...',I7)
      STOP 2
C
C ********** ZERO VARIANCE ABORT
C
  400 WRITE(6,401)
  401 FORMAT('0FUNCTION HAS ZERO VARIANCE')
      STOP 3
   27 FORMAT('0'/'0********'/
     +,' THE DESIRED ACCURACY ',A,' OBTAINED AFTER',I4,' ITERATIONS'//
     + ' FINAL VALUE...',E13.5/
     + ' FINAL ERROR...',E13.5/'0........')
      END
