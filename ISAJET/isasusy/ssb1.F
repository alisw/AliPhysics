#include "isajet/pilot.h"
      COMPLEX*16 FUNCTION SSB1(XS,XMI,XMJ)
C          Modified by M. Drees 10/26/95
#include "isajet/ssinf.inc"
      REAL XS,XMI,XMJ
      DOUBLE PRECISION S,MI,MJ,A0MI,A0MJ
      COMPLEX*16 SSB0
      S=XS
      MI=XMI
      MJ=XMJ
      IF(S.GT.1.D-4*(MI**2+MJ**2)) THEN
        IF(MI.GE.1.D-10) THEN
          A0MI = MI**2*( 1.D0 - LOG(MI**2) + XLAM )
        ELSE
          A0MI = 0.D0
        ENDIF
        IF(MJ.GE.1.D-10) THEN
          A0MJ = MJ**2*( 1.D0 - LOG(MJ**2) + XLAM )
        ELSE
          A0MJ = 0.D0
        ENDIF
        SSB1 = ( (S+MI**2-MJ**2)*SSB0(XS,XMI,XMJ) + A0MJ - A0MI )/2.D0/S
      ELSE IF(ABS(MI-MJ).GT.1.D-4*MJ) THEN
        IF(MI.GT.1.D-10.AND.MJ.GT.1.D-10) THEN
          SSB1 = -(LOG(MI)*(MI**4-2.*MJ**2*MI**2) + MJ**4*LOG(MJ)
     $    -MJ**4/4.D0-.75*MI**4 + MI**2*MJ**2) / (MI**2-MJ**2)**2
     $    + XLAM/2.D0
        ELSEIF(MI.GT.1.D-10) THEN
          SSB1 = -LOG(MI) + .75 + .5*XLAM
        ELSEIF(MJ.GT.1.D-10) THEN
          SSB1 = -LOG(MJ) + .25 + .5*XLAM
        ENDIF
      ELSE IF(MI.NE.0.D0) THEN
        SSB1 = -LOG(MI) + XLAM/2.D0
      ENDIF
      RETURN
      END
