#include "isajet/pilot.h"
C------------------------------------------------------------------
      SUBROUTINE SUGFRZ(Q,G,G0,IG)
C------------------------------------------------------------------
C
C     Freeze out final soft breaking parameters
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sugpas.inc"
      DIMENSION G(29),G0(29)
      INTEGER IG(29)
      REAL Q,MT
      REAL G,G0,TANB
      INTEGER I
C
      TANB=XTANB
      MT=AMT
      DO 200 I=1,5
        G0(I)=G(I)
200   CONTINUE
      IF (Q.LT.MT.AND.IG(6).EQ.0) THEN
        G3MT=G(3)
        G0(6)=G(6)
        IG(6)=1
      END IF
C          Freeze out running gluino mass at MGL
      DO 210 I=7,12
        IF (Q.LT.ABS(G(I)).AND.IG(I).EQ.0) THEN
          G0(I)=G(I)
          IG(I)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=G(I)
        END IF
210   CONTINUE
C          Freeze out Higgs paremeters at HIGFRZ
      DO 211 I=13,14
        IF (Q.LT.HIGFRZ.AND.IG(I).EQ.0) THEN
          G0(I)=G(I)
          IG(I)=1
          G0(I+12)=G(I+12)
          IG(I+12)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=G(I)
          G0(I+12)=G(I+12)
        END IF
211   CONTINUE
C          Freeze out rest at own masses
      DO 220 I=15,24
C        IF (G(I).LT.0.) THEN
C          G(I)=0.
C          NOGOOD=1
C          GO TO 100
C        END IF
        IF (Q.LT.SQRT(ABS(G(I))).AND.IG(I).EQ.0) THEN
          G0(I)=G(I)
          IG(I)=1
        ELSE IF (IG(I).EQ.0) THEN
          G0(I)=G(I)
        END IF
220   CONTINUE
C          Freeze our N_R parameters at Majorana mass scale
      DO 230 I=27,29
        IF (G(I).NE.0.) G0(I)=G(I)
        IF (Q.LT.AMNRMJ.AND.IG(I).EQ.0.) THEN
          IG(I)=1
        END IF
230   CONTINUE
100   RETURN
      END
