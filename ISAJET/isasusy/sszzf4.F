#include "isajet/pilot.h"
      REAL FUNCTION SSZZF4(E)
C-----------------------------------------------------------------------
C          Z_I -> Z_J +B +BBAR VIA HIGGS
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL E
      DOUBLE PRECISION MZ1,MZ2,SN,EE,MH1,MH2,P,MF,BF,TM
C
      MZ2=TMP(1)
      MZ1=TMP(2)
      MH1=TMP(3)
      SN=TMP(4)
      MH2=TMP(5)
      MF=TMP(6)
      EE=E
C
      P=SQRT(MAX(0.D0,EE**2-MZ1**2))
      TM=1.D0-4*MF**2/(MZ2**2+MZ1**2-2*E*MZ2)
      BF=DSQRT(MAX(0.D0,TM))
      SSZZF4=P*BF*(MZ2**2+MZ1**2-2*MZ2*EE-2*MF**2)*
     $(MZ2*EE+SN*MZ2*MZ1)/(MZ2**2+MZ1**2-2*MZ2*EE-MH1**2)/
     $(MZ2**2+MZ1**2-2*MZ2*EE-MH2**2)
      RETURN
      END
