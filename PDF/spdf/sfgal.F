#include "pdf/pilot.h"
       SUBROUTINE SFGAL(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL)
C
C********************************************************************
C*                                                                  *
C*    Parametrization of parton distribution functions              *
C*    in the photon (LO analysis) - full  solution of AP eq.!       *
C*                                                                  *
C* authors: H.Abramowicz, E. Gurvich and A.Levy  (GAL)              *
C*                 /hep-ph/9711355/                                 *
C*                                                                  *
C* Prepared by:                                                     *
C*             Dr E. Gurvich                                        *
C*             Tel Aviv University                                  *
C*             e-mail: gurvich@post.tau.ac.il                       *
C*                                                                  *
C* Modified by:                                                     *
C*             H. Plothow-Besch/CERN-PPE                            *
C*                                                                  *
C********************************************************************
C
C...LAMBDA=0.2 GeV, Nf=4  (i.e. xb=0)
C...range of validity:
C...    10**-4 < X  < 1
C...     4(1)  < Q2 < 10**5 GeV**2
C
#include "pdf/impdp.inc"
#include "pdf/expdp.inc"
     +        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL,
     +        XPDF(-6:6), X, Q, Q2
C
       X = DX
       Q = DQ
       Q2 = Q*Q
       CALL PHGAL(X,Q2,XPDF)
C
      DDV = XPDF(1)
      DUV = XPDF(2)
      DDB = XPDF(-1)
      DUB = XPDF(-2)
      DSB = XPDF(-3)
      DCB = XPDF(-4)
      DBB = XPDF(-5)
      DTB = XPDF(-6)
      DGL = XPDF(0)
C
      RETURN
      END
