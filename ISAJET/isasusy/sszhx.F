#include "isajet/pilot.h"
      FUNCTION SSZHX(X)
C----------------------------------------------------------------
C          Auxiliary function for Z -> HL Z*. Called by SSTEST.
C----------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
#include "isajet/sstype.inc"
C
      REAL X,SSZHX
      REAL MHL,MZ,GZ,R,DEN,TERM
C
      MHL=TMP(1)
      MZ=AMZ
      GZ=GAMZ
      R=MHL/MZ
      TERM=(1.-X+X**2/12.+2*R**2/3.)*SQRT(X**2-4*R**2)
      DEN=(X-R**2)**2+(GZ/MZ)**2
      SSZHX=TERM/DEN
      RETURN
      END
