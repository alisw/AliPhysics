#include "isajet/pilot.h"
      DOUBLE PRECISION FUNCTION SSMQCD(DM,DQ)
C-----------------------------------------------------------------------
C     Calculate leading-log running mass for quark with mass DM at 
C     scale Q, using alpha_s which is continuous across thresholds.
C     See Drees and Hikasa, Phys. Lett. B240: 455-464, Eq. 4.5.
C
C     Note the threshold is at Q = 2 m, not at Q = m as in MSbar.
C
C     Bisset's QCDRAD, WDHFFC
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
C
      DOUBLE PRECISION DM,DQ,DLAM4,DLAM5,DLAM6,DNEFF,POW,RENORM
     $,DQBT,DQTP
C
C         Do nothing for light quarks
C
      IF(DM.LT.1.0) THEN
        SSMQCD=DM
        RETURN
      ENDIF
C
C          Calculate running mass
C
      DLAM4=DBLE(ALQCD4)
      DQBT=2*DBLE(AMBT)
      DQTP=2*DBLE(AMTP)
      SSMQCD=0
C          Q < 2 m(b)
      DNEFF=4
      POW=12.D0/(33.D0-2.*DNEFF)
      IF(DQ.LT.DQBT) THEN
        RENORM=(LOG(2*DM/DLAM4)/LOG(DQ/DLAM4))**POW
        SSMQCD=RENORM*DM
        RETURN
      ELSE
        RENORM=(LOG(2*DM/DLAM4)/LOG(DQBT/DLAM4))**POW
      ENDIF
C          2 m(b) < Q < 2 m(t)
      DNEFF=5
      POW=12.D0/(33.D0-2.*DNEFF)
      DLAM5=DEXP((25.D0*LOG(DLAM4)-LOG(DQBT**2))/23.D0)
      IF(DQ.GE.DQBT.AND.DQ.LT.DQTP) THEN
        RENORM=RENORM
     $  *(LOG(DQBT/DLAM5)/LOG(DQ/DLAM5))**POW
        SSMQCD=RENORM*DM
        RETURN
      ELSE
        RENORM=RENORM
     $  *(LOG(DQBT/DLAM5)/LOG(DQTP/DLAM5))**POW
      ENDIF
C          2 m(t) < Q
      DNEFF=6
      POW=12.D0/(33.D0-2.*DNEFF)
      DLAM6=DEXP((25.D0*LOG(DLAM4)-LOG(DQBT**2)
     $-LOG(4*AMTP**2))/21.D0) 
      RENORM=RENORM
     $*(LOG(DQTP/DLAM6)/LOG(DQ/DLAM6))**POW
      SSMQCD=RENORM*DM
      RETURN
      END
