#include "isajet/pilot.h"
      REAL FUNCTION FBRBM(X)
C
C     Integrand for convolution of 
C     bremsstrahlung with beamstrahlung spectra
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/brembm.inc"
C
      REAL EBEAM,ESTRUC,X
C
      FBRBM=EBEAM(X,EB)*ESTRUC(XMIN/X,QSQBM)/X
      RETURN
      END
