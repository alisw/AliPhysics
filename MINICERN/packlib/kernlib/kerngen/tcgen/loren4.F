*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE LOREN4  (DIR,P4IN,P4OUT)
C
C CERN PROGLIB# U101    LOREN4          .VERSION KERNFOR  4.07  830624
C ORIG. 20/02/68
C
#if !defined(CERNLIB_B48M)
      DOUBLE PRECISION PCM2, ONMCM, EPBETA, PROD
#endif
      DIMENSION DIR(4),P4IN(4),P4OUT(4)
C
C--                VN(A) MEANS N-VECTOR A
C--                GAMMA=ECM/MCM
C--                EPBETA=ECM*V3(PCM)*V3(BETA)
C--                V3(BETA)=V3(PCM)/ECM
C
      PCM2=DIR(1)*DIR(1)+DIR(2)*DIR(2)+DIR(3)*DIR(3)
#if defined(CERNLIB_B48M)
      ONMCM=1.D0/ SQRT (DIR(4)*DIR(4)-PCM2)
#endif
#if !defined(CERNLIB_B48M)
      ONMCM=1.D0/DSQRT (DIR(4)*DIR(4)-PCM2)
#endif
      EPBETA=P4IN(1)*DIR(1)+P4IN(2)*DIR(2)+P4IN(3)*DIR(3)
      PROD=EPBETA*(DIR(4)*ONMCM-1.D0)/PCM2-P4IN(4)*ONMCM
      P4OUT(4)=ONMCM*(P4IN(4)*DIR(4)-EPBETA)
         DO 50 I=1,3
   50 P4OUT(I)=P4IN(I)+DIR(I)*PROD
      RETURN
      END
