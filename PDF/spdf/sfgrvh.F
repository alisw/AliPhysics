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
C
      SUBROUTINE SFGRVH(DX,DSCAL,DUPV,DDNV,DSEA,DSTR,DCHM,DBOT,DTOP,DGL)
C
C :::::::::::: Structure functions from Glueck, Reya and Vogt  (HO)
C :::::::::::: Lambda(NF=3) = 0.200 GeV
C
      DIMENSION XQ(7)
#include "pdf/expdp.inc"
     +       DX,DSCAL,DUPV,DDNV,DSEA,DSTR,DCHM,DBOT,DTOP,DGL
      PARAMETER (ALAM=0.200)
C
      X = DX
      Q2 = DSCAL**2
C
C     generate the individual structure fcn calls
C
C...Higher order parametrization of Gluck, Reya, Vogt.
C...Allowed variable range: 0.2 GeV^2 < Q2 < 1E6 GeV^2; 1E-4 < x < 1.

C...Determine s and x.
          Q2IN=MIN(1.E6,MAX(0.2,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(0.2/ALAM**2))
        SD2=SD**2
        XC=MAX(1.E-4,X)
        XL=-LOG(XC)

C...Calculate structure functions.
        XQ(1)=(1.364+0.989*SD-0.236*SD2)*XC**(0.593-0.048*SD)*
     &  (1.+(8.912-6.092*SD+0.852*SD2)*XC+(-16.737+7.039*SD)*XC**2+
     &  (10.275+0.806*SD-2.000*SD2)*XC**3)*
     &  (1.-XC)**(2.043+1.408*SD-0.283*SD2)
        XQ(2)=(0.835+0.527*SD-0.144*SD2)*XC**(0.600-0.054*SD)*
     &  (1.+(10.245-7.821*SD+1.325*SD2)*XC+(-19.511+10.940*SD-
     &  1.133*SD2)*XC**2+(12.836-2.570*SD-1.041*SD2)*XC**3)*
     &  (1.-XC)**(3.083+1.382*SD-0.276*SD2)
        XQ(3)=(XC**(0.321-0.135*SD)*((10.51-2.299*SD)+
     &  (-17.28+0.755*SD)*XC+(8.242+2.543*SD)*XC**2)*
     &  XL**(-2.023-0.103*SD)+SD**1.044*
     &  EXP(-(-1.178+2.792*SD)+SQRT(2.318*SD**1.673*XL)))*
     &  (1.-XC)**(3.720+2.337*SD-0.199*SD2)
        XQ(4)=SD**0.761*(1.+(6.078-2.065*SD)*XC)*(1.-XC)**(4.654+
     &  0.603*SD-0.326*SD2)*EXP(-(4.231+1.036*SD)+SQRT(3.419*SD**0.316*
     &  XL))/XL**(0.897-0.618*SD)
        XQ(5)=0.
        IF(SD.GT.0.918) XQ(5)=(SD-0.918)*(1.-XC)**(3.328+0.859*SD)*
     &  EXP(-(3.837+1.504*SD)+SQRT((2.150+1.291*SD)*SD**0.5*XL))
        XQ(6)=0.
        IF(SD.GT.1.353) XQ(6)=(SD-1.353)*(1.-XC)**(3.382+0.909*SD)*
     &  EXP(-(4.130+1.486*SD)+SQRT((2.895+1.240*SD)*SD**0.5*XL))

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
