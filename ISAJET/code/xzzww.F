#include "isajet/pilot.h"
      SUBROUTINE XZZWW
C
C          SET UP Z0 Z0 -> W+ W- AMPLITUDES AS RATIONAL FUNCTIONS OF Z
C
C          RE(F(Z,L)) = SUM(I,J)(ANWWWW(I+1,J,L)*Z**I
C                                  /(ADWWWW(1,J)+ADWWWW(2,J)*Z))
C          IM(F(Z,L)) = AIWWWW(L)   (INDEPENDENT OF Z)
C          J LABELS PIECES WITH SAME DENOMINATOR.
C          L=1 FOR 0,0; L=2 FOR 1,-1; L=3 FOR 1,1; L=4 FOR 0,1
C
C          *NOTE* A FACTOR OF SIN(THETA)/SQRT(2) IS REMOVED FROM F01
C
#include "isajet/const.inc"
#include "isajet/jetpar.inc"
#include "isajet/wcon.inc"
#include "isajet/hcon.inc"
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION WM,ZM,ZM2,ZM3,ZM4,ZM5,ZM6,HM,HM2,HM3,HM4,HG,HG2
     $,PROPH,RTS,S,S2,S3,SW,QQ0,QQI,QQF
#endif
C
C          USE UNITS OF WM TO AVOID LARGE NUMBERS - NOTE ANWWWW/ADWWWW
C          AND AIWWWW ARE DIMENSIONLESS
      WM=WMASS(2)
      ZM=WMASS(4)/WM
      ZM2=ZM**2
      ZM3=ZM**3
      ZM4=ZM**4
      ZM5=ZM**5
      ZM6=ZM**6
      HM=HMASS/WM
      HM2=HM**2
      HM3=HM**3
      HM4=HM**4
      HG=HGAM/WM
      HG2=HG**2
      RTS=QMW/WM
      S=RTS**2
      S2=S**2
      S3=S**3
      PROPH=(S-HM2)**2+(HM*HG)**2
C          CORRECT SIGN OF HIGGS AMPLITUDE.
      PROPH=-PROPH
C
      CW=1./ZM
      CW2=CW**2
      SW2=1.-CW2
      SW=SQRT(SW2)
      QQ0=.5*RTS
      QQI=.5*SQRT(S-4.*ZM2)
      QQF=.5*SQRT(S-4.)
      GSQ=4.*PI*ALFA/SW2
C
C          FROM ZZWW3.EX
      ANWWWW(1,1,1) = -1.60E+01 * ((HM2 * ZM3) / (CW * PROPH))
     $ + 1.60E+01 * ((S * ZM3) / (CW * PROPH) ) - 8.00E+00 * ((S2
     $ * ZM) / (CW * PROPH)) - 8.00E+00 * ((S2 * ZM3) / (CW * PROPH))
     $ + 4.00E+00 * ((S3 * ZM) / (CW * PROPH)) + 8.00E+00 * ((HM2 * S
     $ * ZM) / (CW * PROPH)) + 8.00E+00 * ((HM2 * S * ZM3) / (CW
     $ * PROPH)) - 4.00E+00 * ((HM2 * S2 * ZM) / (CW * PROPH))
     $ + 8.00E+00 * CW2 * S - 6.00E+00 * CW2 * S2 + 8.00E+00 * CW2
     $ * S * ZM2
      ANWWWW(1,1,2) = 4.00E+00 * CW2 * S
      ANWWWW(1,1,3) = -1.60E+01 * ((HM2 * ZM3) / (CW * PROPH))
     $ + 1.60E+01 * ((S * ZM3) / (CW * PROPH) ) - 8.00E+00 * ((S2
     $ * ZM) / (CW * PROPH)) + 8.00E+00 * ((HM2 * S * ZM) / (CW
     $ * PROPH)) + 1.20E+01 * CW2 * S - 3.20E+01 * CW2 * ZM2
      ANWWWW(1,1,4) = 0.00E+00
      ANWWWW(1,2,1) = -4.00E+00 * CW2 * S2 + 3.00E+00 * CW2 * S3
     $ - 9.60E+01 * CW2 * ZM4 + 3.20E+01 * CW2 * ZM6 + 8.00E+00 * CW2
     $ * S * ZM2 + 1.60E+01 * CW2 * S * ZM4 - 8.00E+00 * CW2 * S
     $ * ZM6 - 1.00E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(1,2,2) = -4.00E+00 * CW2 * S - 2.00E+00 * CW2 * S2
     $ + 6.40E+01 * CW2 * ZM4 + 1.20E+01 * CW2 * S * ZM2
      ANWWWW(1,2,3) = 4.00E+00 * CW2 * S - 6.00E+00 * CW2 * S2
     $ - 6.40E+01 * CW2 * ZM4 + 2.00E+01 * CW2 * S * ZM2
      ANWWWW(1,2,4) = 6.40E+01 * CW2 * QQ0 * QQF * QQI + 1.60E+01
     $ * CW2 * QQ0 * QQF * QQI * S + 9.60E+01 * CW2 * QQ0 * QQF * QQI
     $ * ZM2 + 3.20E+01 * CW2 * QQ0 * QQF * QQI * ZM4
      ANWWWW(1,3,1) = -4.00E+00 * CW2 * S2 + 3.00E+00 * CW2 * S3
     $ - 9.60E+01 * CW2 * ZM4 + 3.20E+01 * CW2 * ZM6 + 8.00E+00 * CW2
     $ * S * ZM2 + 1.60E+01 * CW2 * S * ZM4 - 8.00E+00 * CW2 * S
     $ * ZM6 - 1.00E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(1,3,2) = -4.00E+00 * CW2 * S - 2.00E+00 * CW2 * S2
     $ + 6.40E+01 * CW2 * ZM4 + 1.20E+01 * CW2 * S * ZM2
      ANWWWW(1,3,3) = 4.00E+00 * CW2 * S - 6.00E+00 * CW2 * S2
     $ - 6.40E+01 * CW2 * ZM4 + 2.00E+01 * CW2 * S * ZM2
      ANWWWW(1,3,4) = -6.40E+01 * CW2 * QQ0 * QQF * QQI - 1.60E+01
     $ * CW2 * QQ0 * QQF * QQI * S - 9.60E+01 * CW2 * QQ0 * QQF * QQI
     $ * ZM2 - 3.20E+01 * CW2 * QQ0 * QQF * QQI * ZM4
      ANWWWW(1,4,1) = 0.00E+00
      ANWWWW(1,4,2) = 0.00E+00
      ANWWWW(1,4,3) = 0.00E+00
      ANWWWW(1,4,4) = 0.00E+00
      ANWWWW(2,1,1) = 0.00E+00
      ANWWWW(2,1,2) = 0.00E+00
      ANWWWW(2,1,3) = 0.00E+00
      ANWWWW(2,1,4) = 8.00E+00 * CW2 * QQ0 * S
      ANWWWW(2,2,1) = 4.80E+01 * CW2 * QQF * QQI * S - 2.00E+01 * CW2
     $ * QQF * QQI * S2 + 6.40E+01 * CW2 * QQF * QQI * ZM2 + 9.60E+01
     $ * CW2 * QQF * QQI * S * ZM2 - 1.60E+01 * CW2 * QQF * QQI * S
     $ * ZM4
      ANWWWW(2,2,2) = 8.00E+00 * CW2 * QQF * QQI * S
      ANWWWW(2,2,3) = 5.60E+01 * CW2 * QQF * QQI * S
      ANWWWW(2,2,4) = -8.00E+00 * CW2 * QQ0 * S2 + 1.92E+02 * CW2
     $ * QQ0 * ZM2 + 1.28E+02 * CW2 * QQ0 * ZM4 - 8.00E+00 * CW2
     $ * QQ0 * S * ZM2 - 8.00E+00 * CW2 * QQ0 * S * ZM4
      ANWWWW(2,3,1) = -4.80E+01 * CW2 * QQF * QQI * S + 2.00E+01 * CW2
     $ * QQF * QQI * S2 - 6.40E+01 * CW2 * QQF * QQI * ZM2
     $ - 9.60E+01 * CW2 * QQF * QQI * S * ZM2 + 1.60E+01 * CW2 * QQF
     $ * QQI * S * ZM4
      ANWWWW(2,3,2) = -8.00E+00 * CW2 * QQF * QQI * S
      ANWWWW(2,3,3) = -5.60E+01 * CW2 * QQF * QQI * S
      ANWWWW(2,3,4) = -8.00E+00 * CW2 * QQ0 * S2 + 1.92E+02 * CW2
     $ * QQ0 * ZM2 + 1.28E+02 * CW2 * QQ0 * ZM4 - 8.00E+00 * CW2
     $ * QQ0 * S * ZM2 - 8.00E+00 * CW2 * QQ0 * S * ZM4
      ANWWWW(2,4,1) = 0.00E+00
      ANWWWW(2,4,2) = 0.00E+00
      ANWWWW(2,4,3) = 0.00E+00
      ANWWWW(2,4,4) = 0.00E+00
      ANWWWW(3,1,1) = 2.00E+00 * CW2 * S2
      ANWWWW(3,1,2) = -4.00E+00 * CW2 * S
      ANWWWW(3,1,3) = 4.00E+00 * CW2 * S
      ANWWWW(3,1,4) = 0.00E+00
      ANWWWW(3,2,1) = 3.20E+01 * CW2 * S - 1.60E+01 * CW2 * S2 + CW2
     $ * S3 + 9.60E+01 * CW2 * S * ZM2 + 3.20E+01 * CW2 * S * ZM4
     $ - 2.20E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(3,2,2) = 4.00E+00 * CW2 * S + 2.00E+00 * CW2 * S2
     $ - 6.40E+01 * CW2 * ZM4 - 1.20E+01 * CW2 * S * ZM2
      ANWWWW(3,2,3) = 2.80E+01 * CW2 * S - 1.00E+01 * CW2 * S2
     $ + 6.40E+01 * CW2 * ZM4 + 1.20E+01 * CW2 * S * ZM2
      ANWWWW(3,2,4) = 1.60E+01 * CW2 * QQ0 * QQF * QQI * S
      ANWWWW(3,3,1) = 3.20E+01 * CW2 * S - 1.60E+01 * CW2 * S2 + CW2
     $ * S3 + 9.60E+01 * CW2 * S * ZM2 + 3.20E+01 * CW2 * S * ZM4
     $ - 2.20E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(3,3,2) = 4.00E+00 * CW2 * S + 2.00E+00 * CW2 * S2
     $ - 6.40E+01 * CW2 * ZM4 - 1.20E+01 * CW2 * S * ZM2
      ANWWWW(3,3,3) = 2.80E+01 * CW2 * S - 1.00E+01 * CW2 * S2
     $ + 6.40E+01 * CW2 * ZM4 + 1.20E+01 * CW2 * S * ZM2
      ANWWWW(3,3,4) = -1.60E+01 * CW2 * QQ0 * QQF * QQI * S
      ANWWWW(3,4,1) = 0.00E+00
      ANWWWW(3,4,2) = 0.00E+00
      ANWWWW(3,4,3) = 0.00E+00
      ANWWWW(3,4,4) = 0.00E+00
      ANWWWW(4,1,1) = 0.00E+00
      ANWWWW(4,1,2) = 0.00E+00
      ANWWWW(4,1,3) = 0.00E+00
      ANWWWW(4,1,4) = 0.00E+00
      ANWWWW(4,2,1) = 4.00E+00 * CW2 * QQF * QQI * S2
      ANWWWW(4,2,2) = -8.00E+00 * CW2 * QQF * QQI * S
      ANWWWW(4,2,3) = 8.00E+00 * CW2 * QQF * QQI * S
      ANWWWW(4,2,4) = 0.00E+00
      ANWWWW(4,3,1) = -4.00E+00 * CW2 * QQF * QQI * S2
      ANWWWW(4,3,2) = 8.00E+00 * CW2 * QQF * QQI * S
      ANWWWW(4,3,3) = -8.00E+00 * CW2 * QQF * QQI * S
      ANWWWW(4,3,4) = 0.00E+00
      ANWWWW(4,4,1) = 0.00E+00
      ANWWWW(4,4,2) = 0.00E+00
      ANWWWW(4,4,3) = 0.00E+00
      ANWWWW(4,4,4) = 0.00E+00
C
      ADWWWW(1,1) = 1.00E+00
      ADWWWW(1,2) = S - 2.00E+00 * ZM2
      ADWWWW(1,3) = S - 2.00E+00 * ZM2
      ADWWWW(1,4) = 1.00E+00
      ADWWWW(2,1) = 0.00E+00
      ADWWWW(2,2) = -1.00E+00 * S + 4.00E+00 * (S / (S + 4.00E+00
     $ * QQF * QQI)) - 1.60E+01 * (ZM2 / (S + 4.00E+00 * QQF * QQI))
     $ + 4.00E+00 * ((S * ZM2) / (S + 4.00E+00 * QQF * QQI))
      ADWWWW(2,3) = S - 4.00E+00 * (S / (S + 4.00E+00 * QQF * QQI))
     $ + 1.60E+01 * (ZM2 / (S + 4.00E+00 * QQF * QQI)) - 4.00E+00
     $ * ((S * ZM2) / (S + 4.00E+00 * QQF * QQI))
      ADWWWW(2,4) = 0.00E+00
C
      AIWWWW(1) = 1.60E+01 * ((HG * HM * ZM3) / (CW * PROPH))
     $ - 8.00E+00 * ((HG * HM * S * ZM) / (CW * PROPH))
     $ - 8.00E+00 * ((HG * HM * S * ZM3) / (CW * PROPH)) + 4.00E+00
     $ * ((HG * HM * S2 * ZM) / (CW * PROPH))
      AIWWWW(2) = 0.00E+00
      AIWWWW(3) = 1.60E+01 * ((HG * HM * ZM3) / (CW * PROPH))
     $ - 8.00E+00 * ((HG * HM * S * ZM) / (CW * PROPH))
      AIWWWW(4) = 0.00E+00
C
C          RESTORE MISSING FACTORS
      DO 100 J=1,4
      AIWWWW(J)=AIWWWW(J)*GSQ/(16.*ZM2)
      DO 100 I=1,4
      DO 110 K=1,4
110   ANWWWW(K,I,J)=ANWWWW(K,I,J)*GSQ/(16.*ZM2)
100   CONTINUE
C
      RETURN
      END
