*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:58  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:07  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C----------------------------------------------------------------------
      FUNCTION VALDIS(XX)
      COMMON/W5051I2/ Q2L
      PARAMETER (Q02=0.5,ALAM2=0.045*0.045)
      PARAMETER (AL0=0.7, BET0=3.21, ALPRM=-0.19, BETPRM=0.82)

      DATA ZERO/0./,HALF/0.5/,ONE/1./,TWO/2./,THRE/3./
      DATA IRAZ/0/, Q2OLD/-1./
      SAVE AV, ALPHA, BETA

      IF(Q2L.NE.Q2OLD) THEN
          Q2OLD = Q2L
          S     = LOG( LOG( MAX( Q02,Q2L)/ALAM2)/ LOG( Q02/ALAM2))
          ALPHA = AL0  + ALPRM *S
          BETA  = BET0 + BETPRM*S
          AV    = THRE*GAMMA(ALPHA+BETA+ONE)/
     1                                 (GAMMA(ALPHA)*GAMMA(BETA+ONE))
      ENDIF
      VALDIS=AV * XX**ALPHA * (ONE-XX)**BETA
      RETURN
      END
