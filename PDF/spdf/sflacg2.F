*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:58  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:38  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       SUBROUTINE SFLACG2(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL)
C
C********************************************************************
C*                                                                  *
C*    Parametrization of parton distribution functions              *
C*    in the photon (LO analysis) - full  solution of AP eq.!       *
C*                                                                  *
C* authors: H.Abramowicz, K.Charchula and A.Levy  (LAC)             *
C*                 /Phys. Lett. B269 (1991) 458/                    *
C*                                                                  *
C* Prepared by:                                                     *
C*             Krzysztof Charchula, DESY                            *
C*             bitnet: F1PCHA@DHHDESY3                              *
C*             decnet: 13313::CHARCHULA                             *
C*                                                                  *
C* Modified by:                                                     *
C*             H. Plothow-Besch/CERN-PPE                            *
C*                                                                  *
C********************************************************************
C
C...ISET   Q20(GEV2)  comments
C...------------------------------------------------------------
C... 2        4        xG(x,Q20)-->const (fixed) at low x
C...------------------------------------------------------------
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
       CALL PHLAC2(X,Q2,XPDF)
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
