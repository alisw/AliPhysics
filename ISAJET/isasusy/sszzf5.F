#include "isajet/pilot.h"
        REAL FUNCTION SSZZF5(SP)
C-----------------------------------------------------------------------
C          SSWZBF: ziss -> zjss f fbar
C          Drees' I26 integrand for higgs-sfermion interference
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL SP
      DOUBLE PRECISION S,M2,M1,MS,E,QS,MUS,BKT
      DOUBLE PRECISION THZ,TERM,MH,DFI,Q,MF,QP
C
      M2=TMP(1)
      M1=TMP(2)
      MS=TMP(3)
      MH=TMP(5)
      THZ=TMP(4)
      MF=TMP(6)
C
      S=SP
      E=(S+M2**2-M1**2)/2.D0/M2
      QS=E**2-S
      Q=DSQRT(MAX(0.D0,QS))
      QP=Q*DSQRT(MAX(0.D0,1.D0-4*MF**2/S))
      MUS=MS**2+S-M1**2-MF**2
      BKT=S*MS**2-MF**2*(M1**2+M2**2)+THZ*M1*M2*(S-2*MF**2)
      TERM=DLOG((M2*(E+QP)-MUS)/(M2*(E-QP)-MUS))
      DFI=(S*QP/2.D0+BKT*TERM/4.D0/M2)/(S-MH**2)
      SSZZF5=DFI
      RETURN
      END
