#include "pdf/pilot.h"
       SUBROUTINE SASDIR(X,Q2,P2,Q02,XPGA)
C...Purpose: to evaluate the direct contribution, i.e. the C^gamma term,
C...as needed in MSbar parametrizations.
      DIMENSION XPGA(-6:6)
      DATA PMC/1.3/, PMB/4.6/, AEM2PI/0.0011614/
 
C...Reset output.
      DO 100 KFL=-6,6
      XPGA(KFL)=0.
  100 CONTINUE
 
C...Evaluate common x-dependent expression.
      XTMP = (X**2+(1.-X)**2) * (-LOG(X)) - 1.
      CGAM = 3.*AEM2PI*X * (XTMP*(1.+P2/(P2+Q02)) + 6.*X*(1.-X))
 
C...d, u, s part by simple charge factor.
      XPGA(1)=(1./9.)*CGAM
      XPGA(2)=(4./9.)*CGAM
      XPGA(3)=(1./9.)*CGAM
 
C...Also fill for antiquarks.
      DO 110 KF=1,5
      XPGA(-KF)=XPGA(KF)
  110 CONTINUE
 
      RETURN
      END
