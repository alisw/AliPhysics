*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RCA (MM, X, IX, Y, IY)
C     ANALYSIS OR SYNTHESIS OF A REAL EVEN FUNCTION.
C
      REAL X(128), Y(128)
      COMMON /D700DT/ N, N2, N4, M, F, RTTWO
      COMMON /FWORK/ W(321)
      DATA MLOC / - 1 /
C
      M = MM + 1
      IF(M .NE. MLOC) CALL D700SU
      MLOC = M
C
      NI = N4
      NO = N4 + N2 + 1
C
      KX = 1
      NFWA = NI
      NLWA = NI + N2
      DO 10 I = NFWA, NLWA
      W(I) = X(KX)
   10 KX = KX + IX
C
      NQ = N2
      DO 170 L = 1, M
      NQ2 = NQ / 2
      NQ2M1 = NQ2 - 1
C
      NO1 = NO
      NO2 = NO + N2
      NI1 = NI
      NI2 = NI + NQ
      W(NO1) = W(NI1) + W(NI2)
      W(NO2) = W(NI1) - W(NI2)
C
      IF(NQ2M1) 50, 40, 20
C
   20 DO 30 IT = 1, NQ2M1
      NOR1 = NO1 + IT
      NOR2 = NO1 + NQ2 + IT
      NIR1 = NI1 + IT
      NIR2 = NI2 - IT
      W(NOR1) = W(NIR1) + W(NIR2)
      W(NOR2) = W(NIR1) - W(NIR2)
   30 CONTINUE
C
   40 NOR1 = NO1 + NQ2
      NIR1 = NI1 + NQ2
      W(NOR1) = W(NIR1) + W(NIR1)
C
   50 NO1 = NO1 + NQ
      NO2 = NO2 - NQ
      NI1 = NI2 + NQ
      NI2 = NI1 + NQ
      IF(NO1 - NO2) 60, 120, 160
C
   60 NC = 0
      NS = N4
      KC = 0
      KS = N4
C
   70 CONTINUE
      W(NO1) = W(NI1) + W(NI2)
      W(NO2) = W(NI1) - W(NI2)
C
      IF(NQ2M1) 110, 100, 80
C
   80 NC = NC + NQ
      NS = NS - NQ
      CC = W(NC)
      SS = W(NS)
      DO 90 IT = 1, NQ2M1
      NOR1 = NO1 + IT
      NOR2 = NO2 + IT
      NIR1 = NI1 + IT
      NIR2 = NI2 - IT
      NOI1 = NOR1 + NQ2
      NOI2 = NOR2 + NQ2
      NII1 = NIR1 + NQ
      NII2 = NIR2 + NQ
      RE = CC * W(NIR2) - SS * W(NII2)
      AI = SS * W(NIR2) + CC * W(NII2)
      W(NOR1) = W(NIR1) + RE
      W(NOR2) = W(NIR1) - RE
      W(NOI1) = + W(NII1) - AI
      W(NOI2) = - W(NII1) - AI
   90 CONTINUE
C
  100 KC = KC + NQ2
      KS = KS - NQ2
      NOR1 = NO1 + NQ2
      NOR2 = NO2 + NQ2
      NIR1 = NI1 + NQ2
      NIR2 = NI2 + NQ2
      W(NOR1) = 2.0 * (W(KC) * W(NIR1) - W(KS) * W(NIR2))
      W(NOR2) = 2.0 * (W(KS) * W(NIR1) + W(KC) * W(NIR2))
C
  110 NO1 = NO1 + NQ
      NO2 = NO2 - NQ
      NI1 = NI2 + NQ
      NI2 = NI1 + NQ
      IF(NO1 - NO2) 70, 120, 160
C
  120 W(NO1) = W(NI1)
C
      IF(NQ2M1) 160, 150, 130
C
  130 DO 140 IT = 1, NQ2M1
      NOR1 = NO1 + IT
      NOI1 = NOR1 + NQ2
      NIR1 = NI + NQ + IT
      NIR2 = NI + NQ + NQ - IT
      W(NOR1) = + W(NIR1)
      W(NOI1) = - W(NIR2)
  140 CONTINUE
C
  150 NOR1 = NO1 + NQ2
      NIR1 = NI + NQ + NQ2
      W(NOR1) = RTTWO * W(NIR1)
C
  160 NT = NI
      NI = NO
      NO = NT
      NQ = NQ2
  170 CONTINUE
C
      KY = 1
      NFWA = NI
      NLWA = NI + N2
      DO 180 I = NFWA, NLWA
      Y(KY) = W(I) * F
  180 KY = KY + IY
C
      RETURN
      END
