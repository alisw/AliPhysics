*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RPA (MM, X, IX, Y, IY)
C     ANALYSIS OF A REAL PERIODIC FUNCTION.
C
      REAL X(128), Y(128)
      COMMON /D700DT/ N, N2, N4, M, F, RTTWO
      COMMON /FWORK/ W(321)
      DATA MLOC / - 1 /
C
      M = MM
      IF(M .NE. MLOC) CALL D700SU
      MLOC = M
C
      NI = N4 - 1
      NO = NI + N
C
      KX = 1
      NFWA = NI + 1
      NLWA = NI + N
      DO 10 I = NFWA, NLWA
      W(I) = X(KX)
   10 KX = KX + IX
C
      NQ = N2
      DO 80 L = 1, M
C
      NO1 = NO
      NO2 = NO + N2
      NI1 = NI
      NI2 = NI + NQ
C
      DO 20 IT = 1, NQ
      NOR1 = NO1 + IT
      NOR2 = NO2 + IT
      NIR1 = NI1 + IT
      NIR2 = NI2 + IT
      W(NOR1) = W(NIR1) + W(NIR2)
      W(NOR2) = W(NIR1) - W(NIR2)
   20 CONTINUE
C
      NC = 0
      NS = N4
      NO1 = NO1 + NQ
      NO2 = NO2 - NQ
      NI1 = NI2 + NQ
      NI2 = NI1 + NQ
      IF(NO1 - NO2) 30, 50, 70
C
   30 NC = NC + NQ
      NS = NS - NQ
      CC = W(NC)
      SS = W(NS)
C
      DO 40 IT = 1, NQ
      NOR1 = NO1 + IT
      NOI1 = NOR1 + N2
      NOR2 = NO2 + IT
      NOI2 = NOR2 + N2
      NIR1 = NI1 + IT
      NII1 = NIR1 + N2
      NIR2 = NI2 + IT
      NII2 = NIR2 + N2
      RE = CC * W(NIR2) - SS * W(NII2)
      AI = SS * W(NIR2) + CC * W(NII2)
      W(NOR1) = W(NIR1) + RE
      W(NOR2) = W(NIR1) - RE
      W(NOI1) = + W(NII1) + AI
      W(NOI2) = - W(NII1) + AI
   40 CONTINUE
C
      NO1 = NO1 + NQ
      NO2 = NO2 - NQ
      NI1 = NI2 + NQ
      NI2 = NI1 + NQ
      IF(NO1 - NO2) 30, 50, 70
C
   50 DO 60 IT = 1, NQ
      NOR1 = NO1 + IT
      NOI1 = NOR1 + N2
      NIR1 = NI1 + IT
      NIR2 = NI2 + IT
      W(NOR1) = W(NIR1)
      W(NOI1) = W(NIR2)
   60 CONTINUE
C
   70 NT = NI
      NI = NO
      NO = NT
      NQ = NQ / 2
C
   80 CONTINUE
C
      KY = 1
      NFWA = NI + 1
      NLWA = NI + N
      DO 90 I = NFWA, NLWA
      Y(KY) = W(I) * F
   90 KY = KY + IY
      Y(KY) = 0.0
C
      RETURN
      END
