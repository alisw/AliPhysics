#include "isajet/pilot.h"
      REAL FUNCTION SSZZF1(X)
C-----------------------------------------------------------------------
C          SSWZBF: ziss -> zjss f fbar
C          Baer's TFUNC
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL X
      DOUBLE PRECISION N,S,Q,D,SQBRKT,TERM1,QS
      DOUBLE PRECISION MZ2,MZ1,M,DSN
C
      MZ2=TMP(1)
      MZ1=TMP(2)
      M=TMP(3)
      DSN=TMP(4)
C
      N=MZ1**2/MZ2**2
      S=M**2/MZ2**2
      Q=X*(1.D0-N)/2.
      QS=Q**2
      D=(S-2*S*Q-N)/((1.D0-2*Q)*(S-2*Q-N))
      IF(D.LE.0.) THEN
        WRITE(LOUT,*) 'ERROR IN SSZZF1: D,S,Q,N=',D,S,Q,N
        SSZZF1=0
        RETURN
      END IF
      SQBRKT=-Q*(1.D0-2*Q-N)/(1.D0-2*Q)-(2*Q-S+N)/2.D0*DLOG(D)
      TERM1=QS*(1.D0-2*Q-N)**2/(1.D0-2*Q-S)**2/(1.D0-2*Q)
     $+DSQRT(N)/2.D0/(1.D0-2*Q-S)*SQBRKT*DSN
      SSZZF1=(1.D0-N)/2.D0*TERM1
      RETURN
      END
