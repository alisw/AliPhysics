#include "isajet/pilot.h"
C-----------------------------------------------------------------
      SUBROUTINE SURG06(T,GY,F)
C-----------------------------------------------------------------
C
C     Right hand side of truncated renormalization group equations
C          dG_i/dT = F_i(G)
C     using a single common scale MSUSY for SUSY mass thresholds.
C     Added GY(7) for neutrino Yukawa coupling on 9/24/99
C     THIS INCLUDES TWO-LOOP YUKAWAS FOR MSSM ONLY 11/18/99
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sugpas.inc"
      REAL T,GY(7),F(7)
      REAL SINB,BETA,COSB,B1,THTOP,B3,B2,TH2LP,B12,B11,B21,B13,
     $MT,TANB,PI,Q,B33,B23,B22,B32,B31
      REAL A11,A12,A13,A21,A22,A23,A31,A32,A33
      REAL C11,C12,C13,C21,C22,C23,C31,C32,C33
      REAL D11,D12,D13,D21,D22,D23,D31,D32,D33
      REAL C1,C2,C3,CP1,CP2,CP3,CPP1,CPP2,CPP3,BSY1,BSY2,BSY3
      INTEGER NSL,NSD,NSH,NSE,NSU,NSQ,NU,NSG,NSW,NH,NN,NE,ND
c      DATA MZ/91.17/
      DATA ND/3/,NE/3/,NN/3/
      DATA B11/7.96/,B12/5.4/,B13/17.6/,B21/1.8/,B22/25./,B23/24./
      DATA B31/2.2/,B32/9./,B33/14./
      DATA A11/5.2/,A12/2.8/,A13/3.6/,A21/6./,A22/6./,A23/2./
      DATA A31/4./,A32/4./,A33/0./
      DATA C11/1.7/,C12/.5/,C13/1.5/,C21/1.5/,C22/1.5/,C23/.5/
      DATA C31/2./,C32/2./,C33/0./
      DATA D11/3.98/,D12/2.7/,D13/8.8/,D21/.9/,D22/5.833/,D23/12./
      DATA D31/1.1/,D32/4.5/,D33/-26./
      DATA C1/.8666667/,C2/3./,C3/5.333333/
      DATA CP1/.4666667/,CP2/3./,CP3/5.333333/
      DATA CPP1/1.8/,CPP2/3./,CPP3/0./
      DATA BSY1/6.6/,BSY2/1./,BSY3/-3./
      TANB=XTANB
      MT=AMT
C-----THESE ARE VALID FROM MZ TO MGUT 
      Q=MGUT*EXP(T)
      PI=4.*ATAN(1.)
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=SQRT(1.-SINB**2)
      IF (Q.GT.MSUSY) THEN
        NSQ=3
        NSU=3
        NSD=3
        NSL=3
        NSE=3
        NSH=2
        NSW=1
        NSG=1
        NH=2
        TH2LP=1.
      ELSE
        NSQ=0
        NSU=0
        NSD=0
        NSL=0
        NSE=0
        NSH=0
        NSW=0
        NSG=0
        NH=1
        TH2LP=0.
      END IF
      IF (Q.GT.MT) THEN
        NU=3
        THTOP=1.
      ELSE
        NU=2
        THTOP=0.
      END IF
      TH2LP=1.
      THTOP=1.
      B1=2.*(17*NU/12.+5*ND/12.+5*NE/4.+NN/4.)/5.+
     $ NSQ/30.+4*NSU/15.+NSD/15.+NSL/10.+NSE/5.+
     $ 1.*NSH/5.+1.*NH/10.
      B2=-22./3.+.5*(NU+ND)+1.*(NE+NN)/6.+
     $ 1.*NSQ/2.+1.*NSL/6.+1.*NSH/3.+1.*NH/6.+4.*NSW/3.
      B3=2.*(NU+ND)/3.+1.*NSQ/3.+1.*NSU/6.+1.*NSD/6.+2.*NSG-11.
      IF (Q.GT.MSUSY) THEN
      F(1)=GY(1)/16./PI**2*(B1*GY(1)**2+TH2LP/16./PI**2*GY(1)**2*
     $(B11*GY(1)**2+B12*GY(2)**2+B13*GY(3)**2-A11*GY(6)**2-A12*GY(5)**2
     $-A13*GY(4)**2))
      F(2)=GY(2)/16./PI**2*(B2*GY(2)**2+TH2LP/16./PI**2*GY(2)**2*
     $(B21*GY(1)**2+B22*GY(2)**2+B23*GY(3)**2-A21*GY(6)**2-A22*GY(5)**2
     $-A23*GY(4)**2))
      F(3)=GY(3)/16./PI**2*(B3*GY(3)**2+TH2LP/16./PI**2*GY(3)**2*
     $(B31*GY(1)**2+B32*GY(2)**2+B33*GY(3)**2-A31*GY(6)**2-A32*GY(5)**2
     $-A33*GY(4)**2))
      ELSE
      F(1)=GY(1)/16./PI**2*(B1*GY(1)**2+TH2LP/16./PI**2*GY(1)**2*
     $(D11*GY(1)**2+D12*GY(2)**2+D13*GY(3)**2-C11*GY(6)**2-C12*GY(5)**2
     $-C13*GY(4)**2))
      F(2)=GY(2)/16./PI**2*(B2*GY(2)**2+TH2LP/16./PI**2*GY(2)**2*
     $(D21*GY(1)**2+D22*GY(2)**2+D23*GY(3)**2-C21*GY(6)**2-C22*GY(5)**2
     $-C23*GY(4)**2))
      F(3)=GY(3)/16./PI**2*(B3*GY(3)**2+TH2LP/16./PI**2*GY(3)**2*
     $(D31*GY(1)**2+D32*GY(2)**2+D33*GY(3)**2-C31*GY(6)**2-C32*GY(5)**2
     $-C33*GY(4)**2))
      ENDIF
      IF (Q.LT.MSUSY) THEN
        F(4)=GY(4)/16./PI**2*(5*GY(4)**2*COSB**2/2.+3*GY(6)**2*SINB**2*
     $   THTOP+3*GY(5)**2*COSB**2-9*GY(1)**2/4.-9*GY(2)**2/4.
     $   -SINB**2*(3*GY(6)**2*THTOP-3*GY(5)**2-GY(4)**2))
        F(5)=GY(5)/16./PI**2*(9*GY(5)**2*COSB**2/2.+3*GY(6)**2*SINB**2*
     $   THTOP/2.+GY(4)**2*COSB**2-GY(1)**2/4.-9*GY(2)**2/4.
     $   -8*GY(3)**2-SINB**2*(3*GY(6)**2*THTOP-3*GY(5)**2-GY(4)**2))
        F(6)=GY(6)/16./PI**2*(9*GY(6)**2*SINB**2/2.*THTOP+
     $   3*GY(5)**2*COSB**2/2.+GY(4)**2*COSB**2-17.*GY(1)**2/20.
     $   -9*GY(2)**2/4.-8*GY(3)**2+COSB**2*
     $   (3*GY(6)**2*THTOP-3*GY(5)**2-GY(4)**2))
      ELSE
        F(4)=GY(4)/16./PI**2*(4*GY(4)**2+3*GY(5)**2+GY(7)**2-
     $        9*GY(1)**2/5.-3*GY(2)**2+
     $       TH2LP/16./PI**2*((CPP1*BSY1+CPP1**2/2.)*GY(1)**4+
     $       (CPP2*BSY2+CPP2**2/2.)*GY(2)**4+(CPP3*BSY3+CPP3**2/2.)*
     $       GY(3)**4+9*GY(1)**2*GY(2)**2/5.+
     $       GY(5)**2*(-.4*GY(1)**2+16*GY(3)**2)+
     $       GY(4)**2*(1.2*GY(1)**2+6*GY(2)**2)-
     $       (3*GY(6)**2*GY(5)**2+9*GY(5)**4+9*GY(5)**2*GY(4)**2+
     $       10*GY(4)**4)))
        F(5)=GY(5)/16./PI**2*(6*GY(5)**2+GY(6)**2*THTOP+
     $        GY(4)**2-7*GY(1)**2/15.-3*GY(2)**2-16*GY(3)**2/3.+
     $       TH2LP/16./PI**2*((CP1*BSY1+CP1**2/2.)*GY(1)**4+
     $       (CP2*BSY2+CP2**2/2.)*GY(2)**4+(CP3*BSY3+CP3**2/2.)*
     $       GY(3)**4+GY(1)**2*GY(2)**2+8*GY(1)**2*GY(3)**2/9.+
     $       8*GY(2)**2*GY(3)**2+.8*GY(6)**2*GY(1)**2+
     $       GY(5)**2*(.4*GY(1)**2+6*GY(2)**2+16*GY(3)**2)+
     $       1.2*GY(4)**2*GY(1)**2-
     $       (22*GY(5)**4+5*GY(6)**2*GY(5)**2+3*GY(5)**2*GY(4)**2+
     $       3*GY(4)**4+5*GY(6)**4)))
        F(6)=GY(6)/16./PI**2*(6*GY(6)**2*THTOP+GY(5)**2+GY(7)**2-
     $        13*GY(1)**2/15.-3*GY(2)**2-16*GY(3)**2/3.+
     $       TH2LP/16./PI**2*((C1*BSY1+C1**2/2.)*GY(1)**4+
     $       (C2*BSY2+C2**2/2.)*GY(2)**4+(C3*BSY3+C3**2/2.)*
     $       GY(3)**4+GY(1)**2*GY(2)**2+136*GY(1)**2*GY(3)**2/45.+
     $       8*GY(2)**2*GY(3)**2+GY(6)**2*(1.2*GY(1)**2+6*GY(2)**2+
     $       16*GY(3)**2)+.4*GY(5)**2*GY(1)**2-
     $       (22*GY(6)**4+5*GY(6)**2*GY(5)**2+5*GY(5)**4+
     $       GY(5)**2*GY(4)**2)))
      END IF
C     THE NEUTRINO YUKAWA IS TO 1-LOOP ONLY FOR NOW...
      IF (Q.GT.AMNRMJ) THEN
        F(7)=GY(7)/16./PI**2*(3*GY(6)**2+GY(4)**2+4*GY(7)**2-
     $        3*GY(2)**2-3*GY(1)**2/5.)
      ELSE
        F(7)=0.
      END IF
      RETURN
      END
