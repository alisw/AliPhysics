*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:19  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:29  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C----------------------------------------------------------------------
      FUNCTION SEADIS(XX)
      DIMENSION
     +       SQ2(2:3)
      COMMON/W5051I1/
     +           GPSI(2:3),GAA(2:3),GPSIA(2:3),GAPSI(2:3),GMINUS(2:3),
     1           AL(2:3),BET(2:3),GPLUS(2:3),
     2           VNQ02(2:3),SNQ02(2:3),GNQ02(2:3)
      COMMON/W5051I2/ Q2L
*- sea quarks' momenta
      PARAMETER (Q02=0.5,ALAM2=0.045*0.045)

      DATA ZERO/0./,HALF/0.5/,ONE/1./,TWO/2./
      DATA PTF/0.25/,PSF/0.75/
      DATA IRAZ/0/, Q2OLD/-1./
      SAVE SQ2, BSQ2

      IF(Q2L.NE.Q2OLD) THEN
          Q2OLD = Q2L
          S     = LOG( LOG( MAX(Q02,Q2L)/ALAM2)/ LOG( Q02/ALAM2))
          DO I=2,3
             D1    = SNQ02(I)*EXP(-GPSI(I)*S)
             D2PLUS=(ONE-AL(I))*(SNQ02(I)+VNQ02(I))-BET(I)*GNQ02(I)
             D2MINU=AL(I)*(SNQ02(I)+VNQ02(I))+BET(I)*GNQ02(I)
             D2    =D2PLUS*EXP(-GPLUS (I)*S)    +
     1              D2MINU*EXP(-GMINUS(I)*S)    -
     2              VNQ02(I)*EXP(-GPSI(I)*S)
             SQ2(I)=PSF*D2   + PTF*D1
          ENDDO
          BSQ2=SQ2(2)/SQ2(3) - TWO
      ENDIF
      SEADIS=SQ2(2) * (ONE-XX)**BSQ2 * (BSQ2+ONE)
      RETURN
      END
