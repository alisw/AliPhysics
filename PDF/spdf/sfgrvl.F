*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:48  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:35  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE SFGRVL(DX,DSCAL,DUPV,DDNV,DSEA,DSTR,DCHM,DBOT,DTOP,DGL)
C
C :::::::::::: Structure functions from Glueck, Reya and Vogt  (LO)
C :::::::::::: Lambda(NF=3) = 0.250 GeV
C
      DIMENSION XQ(7)
#include "pdf/expdp.inc"
     +       DX,DSCAL,DUPV,DDNV,DSEA,DSTR,DCHM,DBOT,DTOP,DGL
      PARAMETER (ALAM=0.250)
C
      X = DX
      Q2 = DSCAL**2
C
C     generate the individual structure fcn calls
C
C...Lowest order parametrization of Gluck, Reya, Vogt.
C...Allowed variable range: 0.2 GeV^2 < Q2 < 1E6 GeV^2; 1E-4 < x < 1;
C...extended to 0.2 GeV^2 < Q2 < 1E8 GeV^2; 1E-6 < x < 1
C...after consultation with the authors.

C...Determine s and x.
          Q2IN=MIN(1E8,MAX(0.2,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(0.2/ALAM**2))
        XC=MAX(1E-6,X)
        XL=-LOG(XC)

C...Calculate structure functions.
        XQ(1)=(0.794+0.312*SD)*XC**(0.427-0.011*SD)*
     &  (1.+(6.887-2.227*SD)*XC+(-11.083+2.136*SD)*XC**2+
     &  (3.900+1.079*SD)*XC**3)*(1.-XC)**(1.037+1.077*SD)
        XQ(2)=(0.486+0.139*SD)*XC**(0.434-0.018*SD)*
     &  (1.+(7.716-2.684*SD)*XC+(-12.768+3.122*SD)*XC**2+
     &  (4.564+0.470*SD)*XC**3)*(1.-XC)**(1.889+1.129*SD)
        XQ(3)=(XC**(0.415+0.186*SD)*((0.786+0.942*SD)+
     &  (5.256-5.810*SD)*XC+(-4.599+5.842*SD)*XC**2)+SD**0.592*
     &  EXP(-(0.398+2.135*SD)+SQRT(3.779*SD**1.250*XL)))*
     &  (1.-XC)**(1.622+1.980*SD)
        XQ(4)=SD**0.448*(1.-XC)**(5.540-0.445*SD)*
     &  EXP(-(4.668+1.230*SD)+SQRT((13.173-1.361*SD)*SD**0.442*XL))/
     &  XL**(3.181-0.862*SD)
        XQ(5)=0.
        IF(SD.GT.1.125) XQ(5)=(SD-1.125)*(1.-XC)**(2.038+1.022*SD)*
     &  EXP(-(4.290+1.569*SD)+SQRT((2.982+1.452*SD)*SD**0.5*XL))
        XQ(6)=0.
        IF(SD.GT.1.603) XQ(6)=(SD-1.603)*(1.-XC)**(2.230+1.052*SD)*
     &  EXP(-(4.566+1.559*SD)+SQRT((4.147+1.268*SD)*SD**0.5*XL))

C...Put into output array - special factor for small x.
        CXS=1.
C       IF(X.LT.1E-4.AND.ABS(XLOWPW-1.).GT.0.01)
C    &  CXS=(1E-4/X)**(XLOWPW-1.)
      UPV = CXS*XQ(1)
      DUPV = UPV
      DNV = CXS*XQ(2)
      DDNV = DNV
      SEA = CXS*XQ(4)
      DSEA = SEA
      STR = CXS*XQ(4)
      DSTR = STR
      CHM = CXS*XQ(5)
      DCHM = CHM
      DBOT = 0.0D0
      DTOP = 0.0D0
      GL = CXS*XQ(3)
      DGL = GL
C
      RETURN
      END
