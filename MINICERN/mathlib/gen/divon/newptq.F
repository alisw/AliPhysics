*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE NEWPTQ (EPS,T,ETA,SFTBND,XLAMDA,U,FU,GU,XMIN,FMIN,XW,FW
     1,XV,FV,A,FA,B,OLDF,B1,SCXBD,E,D,RR,SS,GTEST1,GTEST2,TOL,ILOC,ITEST
     2)
      INTEGER ILOC, ITEST
      DOUBLE PRECISION EPS, T, ETA, SFTBND, XLAMDA, U, FU, GU,
     1   XMIN, FMIN, XW, FW, XV, FV, A, FA, B, OLDF,
     2   B1, SCXBD, E, D, RR, SS, GTEST1, GTEST2, TOL
      DOUBLE PRECISION A1, D1, D2, Q, R, S, T2, XM
      GOTO (10,20,20,230,220),ILOC
 10   ITEST=2
      TOL=T
      T2=TOL+TOL
      IF(U.LE.0.0D+0.OR.XLAMDA.LE.T2.OR.GU.GT.0.0D+0) RETURN
      ITEST=1
      XMIN=0.0D+0
      XW=0.0D+0
      XV=0.0D+0
      A=0.0D+0
      OLDF=FU
      FMIN=FU
      FW=FU
      FV=FU
      FA=FU
      D=U
      SCXBD=EPS*ABS(XLAMDA)+T
      B=XLAMDA+SCXBD
      E=B
      B1=B
      SCXBD=XLAMDA-SCXBD/(1.0D+0+EPS)
      GTEST1=-1.0D-4*GU
      GTEST2=-ETA*GU
      ILOC=2
      GOTO 190
 20   IF(FU.GT.FMIN) GOTO 50
      IF(U.LT.0.0D+0) GOTO 30
      A=0.0D+0
      FA=FMIN
      GOTO 40
 30   B=0.0D+0
 40   XV=XW
      FV=FW
      FW=FMIN
      FMIN=FU
      XMIN=XMIN+U
      A=A-U
      B=B-U
      XV=XV-U
      XW=0.0D+0-U
      SCXBD=SCXBD-U
      TOL=EPS*ABS(XMIN)+T
      T2=TOL+TOL
      GOTO 90
 50   IF(U.GE.0.0D+0) GOTO 60
      A=U
      FA=FU
      GOTO 70
 60   B=U
 70   IF(FU.GT.FW.AND.XW.NE.0.0D+0) GOTO 80
      XV=XW
      FV=FW
      XW=U
      FW=FU
      GOTO 90
 80   IF(FU.GT.FV.AND.XV.NE.0.0D+0.AND.XV.NE.XW) GOTO 90
      XV=U
      FV=FU
 90   XM=5.0D-1*(A+B)
      IF(ABS(XM).LE.T2-5.0D-1*(B-A).OR.XMIN+B.LE.SFTBND.OR.FA-FMIN.LE.
     1ABS(A)*GTEST2.AND.FMIN.LT.OLDF.AND.(ABS(XMIN-XLAMDA).GT.TOL.OR.S
     2CXBD.LT.B)) GOTO 210
      R=0.0D+0
      Q=0.0D+0
      S=0.0D+0
      IF(ABS(E).LE.TOL) GOTO 120
      IF(ILOC.NE.2) GOTO 100
      Q=2.0D+0*(FW-FMIN-XW*GU)
      S=GU*XW*XW
      IF(XMIN.NE.0.0D+0) S=(2.0D+0*(FMIN-FW)+XW*GU)*XW
      GOTO 110
 100  R=XW*(FV-FMIN)
      Q=XV*(FW-FMIN)
      S=R*XW-Q*XV
      Q=2.0D+0*(Q-R)
 110  IF(Q.GT.0.0D+0) S=-S
      IF(Q.LE.0.0D+0) Q=-Q
      R=E
      IF(D.NE.B1.OR.B.LE.SCXBD) E=D
 120  A1=A
      B1=B
      IF(XMIN.NE.A) GOTO 130
      D=XM
      GOTO 160
 130  IF(B.LE.SCXBD) GOTO 140
      D=-4.0D+0*A
      IF(D.GE.B) D=SCXBD
      GOTO 160
 140  D1=A
      D2=B
      IF(ABS(D2).GT.TOL.AND.(XW.LE.0.0D+0.OR.ABS(D1).LE.TOL)) GOTO 1
     150
      U=D1
      D1=D2
      D2=U
 150  U=-D1/D2
      D=5.0D+0*D2*(1.0D-1+1.0D+0/U)/1.1D+1
      IF(U.LT.1.0D+0) D=5.0D-1*D2*SQRT(U)
 160  IF(D.LE.0.0D+0) A1=D
      IF(D.GT.0.0D+0) B1=D
      IF(ABS(S).GE.ABS(5.0D-1*Q*R).OR.S.LE.Q*A1.OR.S.GE.Q*B1) GOTO 1
     170
      D=S/Q
      IF(D-A.GE.T2.AND.B-D.GE.T2) GOTO 180
      D=TOL
      IF(XM.LE.0.0D+0) D=-TOL
      GOTO 180
 170  E=B
      IF(XM.LE.0.0D+0) E=A
 180  ILOC=3
 190  IF(D.LT.SCXBD) GOTO 200
      D=SCXBD
      SCXBD=SCXBD*(1.0D+0+7.5D-1*EPS)+7.5D-1*TOL
 200  U=D
      IF(ABS(D).LT.TOL.AND.D.LE.0.0D+0) U=-TOL
      IF(ABS(D).LT.TOL.AND.D.GT.0.0D+0) U=TOL
      ITEST=1
      RETURN
 210  RR=XMIN
      SS=5.0D-1
      FU=FMIN
 220  IF(ABS(XMIN-XLAMDA).GE.TOL.OR.XMIN.EQ.T) GOTO 230
      XMIN=XLAMDA
      IF(SCXBD.LE.B) GOTO 230
      U=0.0D+0
      ILOC=4
      ITEST=1
      RETURN
 230  IF(XMIN+B.GT.SFTBND) GOTO 240
      ITEST=4
      RETURN
 240  IF(OLDF-FU.LE.GTEST1*XMIN) GOTO 250
      FMIN=FU
      ITEST=0
      RETURN
 250  IF(XMIN.NE.T) GOTO 260
      ITEST=3
      RETURN
 260  XMIN=RR*SS
      SS=SS*SS
      IF(XMIN.LT.T) XMIN=T
      ITEST=1
      U=0.0D+0
      ILOC=5
      RETURN
      END
