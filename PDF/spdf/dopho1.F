*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:06  cernlib
* Version 7.04
*
* Revision 1.2  1996/10/28 14:59:08  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:17  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       SUBROUTINE DOPHO1(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL)
C********************************************************************
C*                                                                  *
C*    Parametrization of parton distribution functions              *
C*    in the photon (LO analysis) - asymptotic solution of AP eq.!  *
C*                                                                  *
C* authors:  D.Duke and H.Owens (DO)                                *
C*           /Phys.Rev. D26 (1982) 1600/                            *
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
#include "pdf/impdp.inc"
#include "pdf/expdp.inc"
     +        CQ(5),
     +        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL
      PARAMETER (ALPEM=7.29927D-3, PI=3.141592D0)
      PARAMETER (ALAM=0.2D0)
      DATA CQ/0.33333D0,0.66666D0,0.33333D0,0.66666D0,0.33333D0/
C
       Q2 = DQ*DQ
       ALAM2=ALAM**2
       FQ=ALPEM/(2.*PI)*LOG(Q2/ALAM2)
C
C...gluons
       POMG=0.194*(1.-DX)**1.03/(DX**0.97)
       DGL=POMG*FQ
C
C...quarks
        POM1=(1.81-1.67*DX+2.16*DX**2)
        POM2=DX**0.7/(1.-0.4*LOG(1.-DX))
        POM3=38.D-4*(1.-DX)**1.82/(DX**1.18)
          DDB=(CQ(1)**2*POM1*POM2+POM3)*FQ
          DDV=DDB
          DUB=(CQ(2)**2*POM1*POM2+POM3)*FQ
          DUV=DUB
          DSB=(CQ(3)**2*POM1*POM2+POM3)*FQ
          DCB=(CQ(4)**2*POM1*POM2+POM3)*FQ
          DBB=(CQ(5)**2*POM1*POM2+POM3)*FQ
       RETURN
       END
